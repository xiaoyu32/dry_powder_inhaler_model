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

#include "fix_charge_gran.h"

#include "fix_scalar_transport_equation.h"
#include "group.h"
#include "stdlib.h"

#include "update.h"
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
#include "math.h"
#include "fix_polarization.h"

#include "efield_model.h"
#include "efield_model_normal.h"
#include "efield_model_screened.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPolarization::FixPolarization(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix efield/polarization needs per particle radius and mass");
  
  
  bool efieldModelSet = false;
  int iarg = 3;
  
  efieldModel = NULL;
  
  while( iarg < narg )
  {   
      if(strcmp(arg[iarg],"model") == 0)
      {
          ++iarg;
          if( iarg >= narg ){
	      fprintf( screen, "%s \n", arg[iarg-1] ); 
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'model'");
	  }
	  
	  ADDEFIELDMODEL(none,EfieldModel,efieldModel)
	  ADDEFIELDMODEL(normal,EfieldModelNormal,efieldModel)
	  ADDEFIELDMODEL(screened,EfieldModelScreened,efieldModel)
	  
	  if( efieldModel ) efieldModelSet = true;
	  else error->fix_error(FLERR,this,"Unknown electric field model. Available models are: 'none', 'normal' and 'screened'!");
	  	  
      } 
      
      ++iarg;
      
  }
  
  
  if( !efieldModelSet ) error->fix_error(FLERR,this,"No electric field model specified: specify with keyword 'model'");

  
  fix_efield = NULL; 
  fix_p = NULL;
  fix_E = NULL;
  fix_Echarge = NULL;
  
  fix_ef_coupling = NULL;
  ef_coupling_flag = 0;
    
}

FixPolarization::~FixPolarization()
{
    if( efieldModel ) delete efieldModel;
}

void FixPolarization::post_create()
{
  
   /*
   fix_efield = static_cast<FixEfieldGran*> ( modify->find_fix_id("wall_gran_efield") );

   if( !fix_efield ){
   
      char* fixarg[3];
   
      fixarg[0] = (char*) "wall_gran_efield";
      fixarg[1] = (char*) "all";
      fixarg[2] = (char*) "efield/gran/particle";

      modify->add_fix(3,fixarg);
   }
   */
  
  fix_p = static_cast<FixPropertyParticle*>(modify->find_fix_property("polarization","property/particle","vector",3,0,this->style,false));
  if(!fix_p)
  {
  
    char* fixarg[11];
    fixarg[0]= (char *)"polarization";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"polarization";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"no";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"no";
    fixarg[8]= (char *)"0.";
    fixarg[9]= (char *)"0.";
    fixarg[10]= (char *)"0.";
    fix_p = modify->add_fix_property_atom(11,fixarg,style);
  } 
   
  fix_E = static_cast<FixPropertyParticle*>(modify->find_fix_property("electricfieldPol","property/particle","vector",3,0,this->style,false));
  if(!fix_E)
  {
    char* fixarg[11];
    fixarg[0]= (char *)"electricfieldPol";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"electricfieldPol";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"no";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"no";
    fixarg[8]= (char *)"0.";
    fixarg[9]= (char *)"0.";
    fixarg[10]= (char *)"0.";
    fix_E = modify->add_fix_property_atom(11,fixarg,style);
  } 
  
  fix_Echarge = static_cast<FixPropertyParticle*>(modify->find_fix_property("electricfieldPol_","property/particle","vector",3,0,this->style,false));
  if(!fix_Echarge)
  {
    char* fixarg[11];
    fixarg[0]= (char *)"electricfieldPol_";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"electricfieldPol_";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"no";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"no";
    fixarg[8]= (char *)"0.";
    fixarg[9]= (char *)"0.";
    fixarg[10]= (char *)"0.";
    fix_Echarge = modify->add_fix_property_atom(11,fixarg,style);
  } 
  

}

void FixPolarization::init()
{
   

   fix_efield = static_cast<FixEfieldGran*> ( modify->find_fix_style("efield/gran",0) );
   
   if (!fix_efield)  error->all(FLERR,"Please create efield/gran!!!!");
   
   //permittivity = (static_cast<FixPropertyGlobal*>(modify->find_fix_property("permittivity","property/global","scalar",0,0,style)))->compute_scalar();
   //susceptibility = (static_cast<FixPropertyGlobal*>(modify->find_fix_property("susceptibility","property/global","scalar",0,0,style)))->compute_scalar();
   
   if(!force->pair_match("gran", 0)) error->all(FLERR,"Please use a granular pair style for fix efield/gran");
   
   efieldModel->connectToProperties( force->registry );
   
   pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
      
   fix_ef_coupling = static_cast<FixPropertyParticle*>(modify->find_fix_property("electricfield","property/particle","vector",3,0,style,false));
   
   if( fix_ef_coupling ) ef_coupling_flag = 1;
   else  ef_coupling_flag = 0;
}

int FixPolarization::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= PRE_FORCE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

// -- sets internal electric field fixes to zero --
void FixPolarization::initial_integrate(int vflag)
{    

    int nlocal = atom->nlocal;
    int * mask = atom->mask;

    for( int i = 0; i < nlocal; ++i )
    {
        
	if( !( mask[i] && groupbit ) ) continue; 
 
       //fprintf( screen, "Polarization = (%e %e %e)...\n", fix_p->array_atom[i][0], fix_p->array_atom[i][1], fix_p->array_atom[i][2] );

        for( int j = 0; j < 3; ++j )
	{
	    fix_E->array_atom[i][j] = 0.0;
	    fix_Echarge->array_atom[i][j] = 0.0;
	}
    }
    
}

void FixPolarization::pre_force(int vflag)
{
    update_echarge();
    update_epolarization();   
}

void FixPolarization::end_of_step()
{

}

double FixPolarization::compute_scalar()
{
    return 0;
}


//compute electric field contribution of the total electric field
void FixPolarization::update_echarge()
{
    
    if( !fix_efield ) return;
    
    int i,j,k,ii,jj,inum,jnum;

    int *ilist,*jlist,*numneigh,**firstneigh;
    int *touch,**firsttouch;
    double *touchfix,*allhist,**firsthist;
    
    double x[3]; 
    double y[3];
    double normal[3];
    double rr;

    //double cof = 1.0/(4*3.1415926*permittivity);
    
    double* charge = fix_efield->fix_charge->vector_atom;
    
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    
    inum = pair_gran->list->inum;
    ilist = pair_gran->list->ilist;
    numneigh = pair_gran->list->numneigh;
    firstneigh = pair_gran->list->firstneigh;
    
    for( int i = 0; i < nlocal; ++i )
    {
        if ( !(mask[i] & groupbit) ) continue;
    
 	for( k = 0; k < 3; ++k )
	{
	    fix_Echarge->array_atom[i][k] = 0.0;
	}     
    }
        
    // loop over neighbors of my atoms
    for(ii = 0; ii < inum; ii++) {
    //for(i = 0; i < nlocal; i++) {    
    	
	i = ilist[ii];
        jlist = firstneigh[i];
	jnum = numneigh[i];
		
	x[0] = atom->x[i][0];
	x[1] = atom->x[i][1];
	x[2] = atom->x[i][2];
	
	for(jj = 0; jj < jnum; jj++) {
	    
	    j = jlist[jj];
      	    j &= NEIGHMASK;

	    if ( !(mask[i] & groupbit) || !(mask[j] & groupbit) ) continue;
	    
	    if( i < j ) continue;
	    
	    y[0] = atom->x[j][0];
	    y[1] = atom->x[j][1];
	    y[2] = atom->x[j][2];
	    
	    rr = 0;
	    
	    for( k = 0;k < 3; ++k )
	    {
		normal[k] = x[k]-y[k];
		rr += normal[k]*normal[k];
	    }
	    
	    rr = sqrt(rr);
	    
	    for( k = 0; k < 3; ++k )
	    {
	        normal[k] /= rr;
	    }
	    
	    
	    //fprintf(screen, "rr is %e   charge[i] = %e \n", rr, charge[i]);
	    
	    efieldModel->computeEfieldCharge( normal, rr, charge[j], fix_Echarge->array_atom[i] );
	    efieldModel->computeEfieldCharge( normal, rr, -charge[i], fix_Echarge->array_atom[j] );
    
    
	    //fprintf( screen, "i = %d    j = %d   E = ( %e %e %e ) \n", i, j , fix_Echarge->array_atom[i][0], fix_Echarge->array_atom[i][1], fix_Echarge->array_atom[i][2] );
    
	}
       
    }
    
    
    
    if( ef_coupling_flag != 0 )
    {
    
	//large scale electric field contribution
	for( i = 0; i < nlocal; ++i )
	{

	    for( k = 0; k < 3; ++k )
	    {
		fix_Echarge->array_atom[i][k] += fix_ef_coupling->array_atom[i][k]; 
	    }
//fprintf( screen, "i = %d   E = ( %e %e %e ) \n", i,  fix_Echarge->array_atom[i][0], fix_Echarge->array_atom[i][1], fix_Echarge->array_atom[i][2] );
    	  

	}
    
    }
    
    fix_Echarge->do_forward_comm();
    
}

void FixPolarization::update_epolarization()
{
    
    int i,j,k,ii,jj,inum,jnum;

    int *ilist,*jlist,*numneigh,**firstneigh;
    int *touch,**firsttouch;
    double *touchfix,*allhist,**firsthist;
    
    double x[3]; 
    double y[3];
    double normal[3];
    double rr;

    //double cof = 1.0/(4*3.1415926*permittivity);
    
    int * mask = atom->mask;
    int nlocal = atom->nlocal;
    
    inum = pair_gran->list->inum;
    ilist = pair_gran->list->ilist;
    numneigh = pair_gran->list->numneigh;
    firstneigh = pair_gran->list->firstneigh;
    
    for( i = 0; i < nlocal; ++i )
    {
        for( k = 0; k < 3; ++k )
	{
	    fix_E->array_atom[i][k] = 0.0;
	}

    }
    
    
    // loop over neighbors of my atoms
    for(ii = 0; ii < inum; ii++) {
    
	i = ilist[ii];
        jlist = firstneigh[i];
        jnum = numneigh[i];	

	x[0] = atom->x[i][0];
	x[1] = atom->x[i][1];
	x[2] = atom->x[i][2];

	for(jj = 0; jj < jnum; jj++) {
	    
	    j = jlist[jj];
      	    j &= NEIGHMASK;
	    	    
	    if ( !(mask[i] & groupbit) && !(mask[j] & groupbit) ) continue;
	    
	    y[0] = atom->x[j][0];
	    y[1] = atom->x[j][1];
	    y[2] = atom->x[j][2];
	    
	    rr = 0;
	    
	    for( k = 0;k < 3; ++k )
	    {
		normal[k] = y[k]-x[k];
		rr += normal[k]*normal[k];
	    }
	    
	    rr = sqrt(rr);
	    
	    for( k = 0; k < 3; ++k )
	        normal[k] /= rr;

	    efieldModel->computeEfieldPolarization( normal, rr, fix_p->array_atom[j], fix_E->array_atom[i] );
	    efieldModel->computeEfieldPolarization( normal, rr, fix_p->array_atom[i], fix_E->array_atom[j] );	
		
	    
	}
       
    }
    
    
    
    fix_E->do_forward_comm();
    
}


void FixPolarization::get_electricfield( int i, double* p )
{
    
    int nlocal = atom->nlocal;
    
    if( i >= nlocal ) return;
    if( i < 0 ) return;
    
    for( int k = 0; k < 3; ++k )
    {
        p[k] = fix_E->array_atom[i][k] + fix_Echarge->array_atom[i][k];
    }
    
}




//integrate the permanent polarization orientation
void FixPolarization::final_integrate()
{}











