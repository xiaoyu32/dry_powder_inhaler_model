/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Evan Smuts (U Cape Town, surface velocity rotation)
------------------------------------------------------------------------- */

#include "fix_mesh_surface_charge.h"
#include <stdio.h>
#include "string.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "comm.h"
#include "math_extra.h"
#include "fix_property_global.h"
#include "fix_gravity.h"
#include "mpi.h"

#include "fix_charge_gran.h"

#include "efield_model.h"
#include "efield_model_normal.h"
#include "efield_model_screened.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMeshSurfaceCharge::FixMeshSurfaceCharge(LAMMPS *lmp, int narg, char **arg)
: FixMeshSurface(lmp, narg, arg),
  surfaceChargeDensity_(NULL),
  initC( 0 ),
  initFlag( false ),
  new_bd_model_flag( false), 
  triMesh_( NULL ),
  chargeTransferActive( false ),
  fix_charge( NULL ),
  fix_efieldFlux( NULL ),
  fix_directionalEfieldFlux( NULL ),
  fix_prev_loc( NULL ),
  fix_tribo( NULL ),
  forceActive( true ),
  fix_p( NULL ),
  fix_nparcel( NULL )
{
    
    // parse further args 
    bool hasargs = true;
    
    while(iarg_ < narg && hasargs)
    {
    
      hasargs = false;
      if (strcmp(arg[iarg_],"initial_surface_density") == 0) 
      {
	  
	  if (narg <= iarg_+1) error->fix_error(FLERR,this,"not enough arguments");
	  iarg_++;
	  
	  initC = atof( arg[iarg_] );
	  initFlag = true;
	  hasargs = true;  
	  
      }else if (strcmp(arg[iarg_],"force") == 0) 
      {
	  
	  if (narg <= iarg_+1) error->fix_error(FLERR,this,"not enough arguments");
	  iarg_++;
	  
	  if( strcmp( arg[iarg_], "yes" ) == 0 )
	  {
	      forceActive = true;
	  }else{
	     
	      forceActive = false;
	      
	      int rank = 0;
              MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	      if( rank == 0 ) fprintf( screen, "Insulating wall electrostatic force disabled! \n" );
	      
	  }
	  	  
	  hasargs = true;  
	  
      }else if (strcmp(arg[iarg_],"bd_model") == 0)
      {
      	 iarg_++;
	 
	 if( iarg_ >= narg ){
	      fprintf( screen, "%s \n", arg[iarg_-1] ); 
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'bd_model'");
	 }
	 
	 int bd_model = atoi(arg[iarg_]);
	 
	 
	 if (bd_model == 1){
	 	new_bd_model_flag = true;
		
		fprintf( screen, "using new breakdown model...");
	 } 
	 else{
	 	new_bd_model_flag = false;
	 }
		
		
      }  
      
    }
    
    // -- pointer to the mesh object --
    triMesh_ = dynamic_cast<TriMesh*>( this->mesh() );
    
    // -- set wall cfd-dem communication on --
    triMesh_->setActive( true );
    
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceCharge::~FixMeshSurfaceCharge()
{}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceCharge::post_create()
{
    FixMeshSurface::post_create(); 
    mesh()->prop().addElementProperty<ScalarContainer<double> >("surfaceChargeDensity","comm_exchange_borders","frame_invariant","restart_yes"); 
    mesh()->prop().addElementProperty<ScalarContainer<double> >("chargeFlux","comm_exchange_borders","frame_invariant","restart_no");     
}

void FixMeshSurfaceCharge::init()
{

    int rank = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    surfaceChargeDensityFlux_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("chargeFlux");
    surfaceChargeDensityFlux_->setAll( 0.0 );
    
    surfaceChargeDensity_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("surfaceChargeDensity");
    if( initFlag ) surfaceChargeDensity_->setAll( initC );
    
    if( rank == 0 ) fprintf( screen, "Initialized surface charge density.\n" );
    
    fix_tribo = static_cast<FixEfieldGran*>( lmp->modify->find_fix_style("efield/gran", 0) );
    
    if( !fix_tribo )
    {
       error->all( FLERR, "Fix type 'mesh/surface/charge' requires fix of type 'efield/gran' to exist!" );
    }
    
    fix_prev_loc = static_cast<FixPropertyParticle*>(modify->find_fix_property("prev_loc","property/particle","vector",0,0,style));
    fix_charge = static_cast<FixPropertyParticle*>(modify->find_fix_property("charge","property/particle","scalar",0,0,style));
    fix_efieldFlux = static_cast<FixPropertyParticle*>(modify->find_fix_property("efieldFlux","property/particle","scalar",0,0,style));
    fix_directionalEfieldFlux = static_cast<FixPropertyParticle*>(modify->find_fix_property("directionalEfieldFlux","property/particle","vector",0,0,style));
    fix_p = static_cast<FixPropertyParticle*>(modify->find_fix_property("polarization","property/particle","vector",3,0,this->style,false));
    
    int max_type = atom->get_properties()->max_type();
    nparcel = new double[max_type];
    fix_nparcel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("nparcel","property/particle","peratomtype",max_type,0,this->style,false));
    
    for( int i = 0; i < max_type; ++i )
    {
	if( !fix_nparcel ) 
    	{
    		nparcel[i] = 1;
    	}
	else
	{
		nparcel[i] = fix_nparcel->compute_vector(i);

	  	if( nparcel[i] <= 0 ) error->all(FLERR,"Number of particles in a parcel has to be a positive number!");
	}
    }
   
}





double FixMeshSurfaceCharge::calcNearestPoint( const double* x, double* xMesh, int iTri ) const
{
    
    // -- triangle corner points --
    double p1[3] = {0,0,0};
    double p2[3] = {0,0,0};
    double p3[3] = {0,0,0};    
    
    triMesh_->node(iTri,0,p1);
    triMesh_->node(iTri,1,p2);    
    triMesh_->node(iTri,2,p3);
    
    double e1[3] = {0,0,0};
    double e2[3] = {0,0,0};
    double e3[3] = {0,0,0};
    
    minus( p1,p2,e1 );
    minus( p1,p3,e2 );
    
    normalize( e1 );
    normalize( e2 );
    
    cross( e1, e2, e3 );
    
    normalize( e3 );
    
    
    double xp[3] = {0,0,0};
    minus( p1, x, xp );
    
    // -- distance to plane spanned by the face --
    double d = dot( xp,e3 );
    
    // -- x projected on the plane --
    for( int i = 0; i < 3; ++i )
       xMesh[i] = x[i] - d*e3[i];
    
    // -- check if the projected point lies on the face --
    
    // -- compute sum of sub triangles with one corner on the projected point --
    // -- if the sum of the sub triangles is the same as the face area, point is on the face --
    double area[3] = {0,0,0};
    double areaSum = 0;
    
    // -- edges going from xMesh (projected point) to each corner of the triangle -- 
    minus( xMesh, p1 , e1 );
    minus( xMesh, p2 , e2 );
    minus( xMesh, p3 , e3 );
    
    // -- first subtriangle --
    cross( e1,e2,area );
    areaSum += mag( area )/2.0;
    
    // -- second subtriangle --
    cross( e2,e3,area );
    areaSum += mag( area )/2.0;
    
    // -- third subtriangle --
    cross( e3,e1,area );
    areaSum += mag( area )/2.0;

    if( areaSum > ( calcArea( iTri ) + 1e-12 ) ) // -- sum is never smaller than the actual area --
    {
        // -- xMesh does not lie inside the triable, and there is no contact --
	d = centerDistance( x, iTri );
    }
    
    return d>0 ? d:-d;
    
}

double FixMeshSurfaceCharge::centerDistance( const double* x, int iTri ) const
{
    double meshCenter[3] = {0,0,0};

    triMesh_->center( iTri, meshCenter );
    double d = distance( meshCenter, x );
    
    return relu( d );
}

void FixMeshSurfaceCharge::normalVector( const double* x, double* n, int iTri ) const
{
    double meshCenter[3] = {0,0,0};

    triMesh_->center( iTri, meshCenter );
    double d = distance( meshCenter, x );
    
    // -- n pointing from face center to particle center --
    for( int i = 0; i < 3; ++i )
        n[i] = ( x[i] - meshCenter[i] )/d;
        
}

double FixMeshSurfaceCharge::calcArea( int ii ) const
{
    
    double area = 0;
    
    // -- triangle corner points --
    double p1[3];
    double p2[3];
    double p3[3];
    
    // -- triangle edge vectors --
    double a[3];
    double b[3];
    
    double c[3];
    
    if( !this->triMesh_ ) return 0;
    
    triMesh_->node(ii,0,p1);
    triMesh_->node(ii,1,p2);
        
    for( int i = 2; i < this->triMesh_->numNodes(); ++i )
    {    
	triMesh_->node(ii,i,p3);
	minus( p2, p1, a );
	minus( p3, p1, b );
	
	cross(a,b,c);
	
	area += mag( c )/2.0;
	vectorCopy3D( p3, p2 );
    }
    
    return area;
    
}

/* ---------------------------------------------------------------------- */

int FixMeshSurfaceCharge::setmask()
{
    int mask = FixMeshSurface::setmask();
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceCharge::setup_pre_force(int vflag)
{
    FixMeshSurface::setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceCharge::pre_force(int vflag)
{
    FixMeshSurface::pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceCharge::final_integrate()
{
    
    //fprintf( screen, "FixMeshSurfaceCharge::final_integrate() \n" );
    
    FixMeshSurface::final_integrate();
    if( !chargeTransferActive ) return;
    
    //fprintf( screen, "FixMeshSurfaceCharge::final_integrate() %d \n", surfaceChargeDensity_->size() );
    
    // -- add surface charge density fluxes to the face elements --
    for( int i = 0; i < surfaceChargeDensity_->size(); ++i )
    {
        double& chargeFlux = surfaceChargeDensityFlux(i);
        if( chargeFlux != 0 )
	{
	   surfaceChargeDensity(i) += chargeFlux;  
	   chargeFlux = 0;
	}  
    }
    
}


/* ----------------------------------------------------------------------
   return total number of API particles or torque component on body
------------------------------------------------------------------------- */

double FixMeshSurfaceCharge::compute_vector( int n )
{
   return 0;
}

// -- wall face electrostatic force contribution to particles --
void FixMeshSurfaceCharge::compute_scalar_transport( LIGGGHTS::ContactModels::CollisionData& cdata, double * v_wall, int iMesh, int iTri )
{
        
    if( !forceActive ) return;
    
    int ipart = cdata.i;
    
    if( !(atom->mask[ipart] & groupbit) ) return;
    
    int *type = atom->type;

    // -- do not calculate charge transfer for ghosts --
    if( ipart >= atom->nlocal ) return;
    
    const double chargePart = (fix_charge->vector_atom)[ipart];
    
    
    const double faceArea = calcArea( iTri );
    double idt = 1.0/update->dt; 
    
    const double permittivity = fix_tribo->permittivity;
    const double delta_charge = fix_tribo->delta_charge;
    const double bd_field     = fix_tribo->bd_field; 
       
    
    const double chargeFace = surfaceChargeDensity( iTri ) * faceArea;
    
    // -- effective radius of sphere that matches the face area -- (used for charge calculations)
    const double faceReff = sqrt( faceArea / (4 * M_PI ) ); 
    
    // -- compute electric field at the particle position due to the face --
    double Efp[3] = {0,0,0};
    
    double gradEfp[3][3] = { 
    			      {0,0,0},
    			      {0,0,0},
			      {0,0,0}
			   };
    
    // -- particle position --
    const double* x = atom->x[ipart];
    
    // -- closest position on face or face center if not touching --
    double xMesh[3] = {0,0,0};
    double n[3] = {0,0,0};
    
    double d = centerDistance( x, iTri );
    
    // -- cap at the particle radius --
    if( d < atom->radius[ipart] ) d = atom->radius[ipart];
    
    // -- at contact simplifies to infinite wall formula --
    const double deltaf = faceReff + d - atom->radius[ipart];
    
    normalVector( x, n,  iTri );
    
    fix_tribo->efieldModel()->computeEfieldCharge( n, deltaf, chargeFace, Efp );
    if( fix_p ) fix_tribo->efieldModel()->computeGradEfieldCharge( n, deltaf, chargeFace, gradEfp );
    
    double tor[3] = {0,0,0};
    
    // -- electrostatic torque working on the dipole --
    if( fix_p )
    {
        cross( fix_p->array_atom[ipart], Efp, tor );
    }
    
    
    double fn[3] = {0,0,0};
    for ( int i = 0; i < 3; ++i)
    {
    	fn[i] = Efp[i] * chargePart;
	
	if (fix_p)
	{
	     for (int j = 0; j < 3; ++j)
	     {
	     	fn[i] += gradEfp[i][j] * fix_p->array_atom[ipart][j];
	     }
	     atom->torque[ipart][i] += tor[i];
	}
	
	// parcel approach
	fn[i] *= nparcel[ipart];
	
	atom->f[ipart][i] += fn[i];
    }
    
    /*
    for( int i = 0; i < 3; ++i )
    {
        // -- Coulombic force --
	atom->f[ipart][i] += Efp[i] * chargePart;
	
	if( fix_p  )
	{
	    
	    // -- dielectrophoretic force --
	    for( int j = 0; j < 3; ++ j )
	    {
		atom->f[ipart][i] += gradEfp[i][j] * fix_p->array_atom[ipart][j];
	    }
	    
	    atom->torque[ipart][i] += tor[i];
	    
	}
	
    }
    */
    
}
















