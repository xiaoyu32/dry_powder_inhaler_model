/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Efield
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "fix_charge_force_parcel.h"
#include "fix_charge_gran.h"
#include "fix_scalar_transport_equation.h"
#include "group.h"
#include "stdlib.h"
#include "math.h"

#include "atom.h"
#include "compute_pair_gran_local.h"
#include "fix_property_particle.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "mech_param_gran.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair_gran.h"
#include "mpi_liggghts.h"

#include "efield_model.h"
#include "efield_model_normal.h"
#include "efield_model_screened.h"
#include "string.h"
using namespace LAMMPS_NS;
using namespace FixConst;

FixEfieldForceParcel::FixEfieldForceParcel(LAMMPS *lmp, int narg, char **arg) : 
Fix(lmp, narg, arg),
efieldModel_( NULL ),
fix_p( NULL ),
fix_charge( NULL ),
no_torque( true ),
nparcel( NULL ),
fix_nparcel( NULL )
{
    if ((!atom->radius_flag)||(!atom->rmass_flag)) 
      	error->all(FLERR,"Fix efield/force needs per particle radius and mass");
	
	
    int iarg = 3;
    
    while( iarg < narg )
    { 
       if(strcmp(arg[iarg],"model") == 0)
       {
           ++iarg;
           if( iarg >= narg ){
	       fprintf( screen, "%s \n", arg[iarg-1] ); 
	       error->fix_error(FLERR,this,"not enough arguments for keyword 'model'");
	   }

	   ADDEFIELDMODEL(none,EfieldModel,efieldModel_)
	   ADDEFIELDMODEL(normal,EfieldModelNormal,efieldModel_)
	   ADDEFIELDMODEL(screened,EfieldModelScreened,efieldModel_)

	   if( !efieldModel_ ) 
	     error->fix_error
	     (
		  FLERR,
		  this,
		  "Unknown electric field model. Available models are: 'none', 'normal' and 'screened'!"
	      );

       }
       
       ++iarg;
           
    }
    
    if( !efieldModel_ ) error->fix_error(FLERR,this,"No electric field model specified: specify with keyword 'model'");
    
    
    		
}

FixEfieldForceParcel::~FixEfieldForceParcel()
{
    if( efieldModel_ ) delete[] efieldModel_;
}


void FixEfieldForceParcel::init()
{
    pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));

    efieldModel_->connectToProperties( force->registry );
    
    int max_type = atom->get_properties()->max_type();
    fix_nparcel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("nparcel","property/particle","peratomtype",max_type,0,"efield/force/parcel"));

    if( !fix_nparcel ) error->all(FLERR,"Fix 'efield/force/parcel' requires property/particle fix 'nparcel'!");

    nparcel = new double[max_type];
    for( int i = 0; i < max_type; ++i )
    {
    	nparcel[i] = fix_nparcel->compute_vector(i);
	if( nparcel[i] <= 0 ) error->all(FLERR,"Number of particles in a parcel has to be a positive number!");
    }
    
    fix_charge = static_cast<FixPropertyParticle*>(modify->find_fix_property("charge","property/particle","scalar",0,0,this->style,false));
    fix_p = static_cast<FixPropertyParticle*>(modify->find_fix_property("polarization","property/particle","vector",3,0,this->style,false));
    
    if( !fix_charge && !fix_p ) error->fix_error(FLERR,this,"Fix 'efield/force' requires either particle charge or electric dipole!"); 

    fix_init_charge = static_cast<FixPropertyGlobal*>(modify->find_fix_property("initial_charge","property/particle","peratomtype",max_type,0,this->style,false));
    
}

int FixEfieldForceParcel::setmask()
{
    int mask = 0;
    mask |= POST_FORCE;
    return mask;
}

// -- force calculations --
void FixEfieldForceParcel::post_force(int vflag)
{

  int newton_pair = force->newton_pair;
  int inum = pair_gran->list->inum;
  int * ilist = pair_gran->list->ilist;
  int * numneigh = pair_gran->list->numneigh;
  int ** firstneigh = pair_gran->list->firstneigh;
  int *type = atom->type;
  
  double *radius = atom->radius;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  
  for(int ii = 0; ii < inum; ii++) 
  {

    int i = ilist[ii];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    double radi = radius[i];
    int jnum = numneigh[i];
    int * jlist = firstneigh[i];
    int itype = type[i];
    
    for(int jj = 0; jj < jnum; jj++) 
    {
	int j = jlist[jj];
	j &= NEIGHMASK;  
	
	if ( !(mask[i] & groupbit) || !(mask[j] & groupbit) ) continue;
	
	if( i < j ) continue;
	
	const double radi = radius[i];
	const double radj = radius[j]; 
	int jtype = type[j];
	
	double n[3] = {0,0,0};
	double r = 0;
	
	for( int k = 0; k < 3; ++k )
	{
	   n[k] = x[j][k]-x[i][k];
	   r += n[k]*n[k];
	}
	
	// -- center to center distance --
	r = sqrt(r);
	
	// -- compute normal vector --
	for( int k = 0; k < 3; ++k )
	{
	   n[k] /= r;
	}
	
	// -- cap at contact --
	if( r < radi+radj ) r = radi+radj;

	// -- compute electric field at the particle j due to particle i --
	double Efp[3] = {0,0,0};
	
	
	// -- compute electric field gradient at the particle j due to particle i --
	double gradEfp[3][3] = { 
    				  {0,0,0},
    				  {0,0,0},
				  {0,0,0}
			       };
	
	// -- electric field due to particle charge --
	if( fix_charge ) efieldModel()->computeEfieldCharge( n, r, fix_charge->vector_atom[i], Efp );
	if( fix_init_charge ) efieldModel()->computeEfieldCharge( n, r, fix_init_charge->compute_vector(itype-1), Efp);

	if( fix_charge && fix_p ) efieldModel()->computeGradEfieldCharge( n, r, fix_charge->vector_atom[i], gradEfp );	
	
	
	// -- electric field due to particle dipole --
	if( fix_p ) efieldModel()->computeEfieldPolarization( n, r, fix_p->array_atom[i], Efp );
	if( fix_p ) efieldModel()->computeGradEfieldPolarization( n, r, fix_p->array_atom[i], gradEfp );  
	
	
	double fn[3] = {0,0,0};
	double tor[3] = {0,0,0};
	
	double qj; 
	if( fix_charge ) qj = fix_charge->vector_atom[j];
	if( fix_init_charge )  qj = fix_init_charge->compute_vector(jtype-1); 
	
	if( !no_torque && fix_p )
	{
            cross( fix_p->array_atom[j], Efp, tor );
	}
	
	for( int k = 0; k < 3; ++k )
	{
	    // -- Coulombic force --
	    fn[k] += Efp[k] * qj;
	    
	    if( fix_p  )
	    {
	    
		// -- dielectrophoretic force --
		for( int kk = 0; kk < 3; ++ kk )
		{
		    fn[k] += gradEfp[k][kk] * fix_p->array_atom[j][kk];
		}

	    }		
	}
	
	// -- parcel contribution --
	for ( int k = 0; k < 3; ++k)
	{
     	    fn[k] *= nparcel[itype-1] > nparcel[jtype-1]? nparcel[itype-1]: nparcel[jtype-1];
	}
	
	// -- add force to the particles -- (fn is force acting on particle j)
	atom->f[i][0] -= fn[0];
	atom->f[i][1] -= fn[1];
	atom->f[i][2] -= fn[2];
	
	if( !no_torque )
	{
	   atom->torque[i][0] -= tor[0];
	   atom->torque[i][1] -= tor[1];
	   atom->torque[i][2] -= tor[2];
	}
	
	if( j < nlocal )
	{

	   // -- add force to the particles -- (fn is force acting on particle j)
	   atom->f[j][0] += fn[0];
	   atom->f[j][1] += fn[1];
	   atom->f[j][2] += fn[2];

	   if( !no_torque )
	   {
	      atom->torque[i][0] -= tor[0];
	      atom->torque[i][1] -= tor[1];
	      atom->torque[i][2] -= tor[2];
	   }

	}
	
	
    }
    
  }
  
  
}
