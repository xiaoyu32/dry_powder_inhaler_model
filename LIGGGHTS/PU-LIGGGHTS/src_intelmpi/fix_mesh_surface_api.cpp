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

#include "fix_mesh_surface_api.h"
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

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMeshSurfaceApi::FixMeshSurfaceApi(LAMMPS *lmp, int narg, char **arg)
: FixMeshSurface(lmp, narg, arg),
  surface_density_(NULL),
  initC( 0 ),
  initFlag( false ),
  triMesh_( NULL )
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
	  
      } 
      
    }
    
    // -- pointer to the mesh object --
    triMesh_ = dynamic_cast<TriMesh*>( this->mesh() );
    
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceApi::~FixMeshSurfaceApi()
{

}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceApi::post_create()
{
    FixMeshSurface::post_create(); 
    mesh()->prop().addElementProperty<ScalarContainer<double> >("apiSurfaceDensity","comm_exchange_borders","frame_invariant","restart_yes");   
}

void FixMeshSurfaceApi::init()
{

    int rank = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    surface_density_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("apiSurfaceDensity");
    if( initFlag ) surface_density_->setAll( initC );
    
    if( rank == 0 ) fprintf( screen, "Initialized API wall concentration values.\n" );
    
    // -- debug --
    surface_density_->setAll( double( rank+1 ) );
    
}



double FixMeshSurfaceApi::calcArea( int ii ) const
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

int FixMeshSurfaceApi::setmask()
{
    int mask = FixMeshSurface::setmask();
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceApi::setup_pre_force(int vflag)
{
    FixMeshSurface::setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceApi::pre_force(int vflag)
{
    FixMeshSurface::pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceApi::final_integrate()
{
    FixMeshSurface::final_integrate();
}


/* ----------------------------------------------------------------------
   return total number of API particles or torque component on body
------------------------------------------------------------------------- */

double FixMeshSurfaceApi::compute_vector(int n)
{
   return 0;
}

void FixMeshSurfaceApi::minus( const double* a, const double* b, double* c ) const
{
   for( int i = 0; i < 3; ++i )
      c[i] = b[i]-a[i];
}


void FixMeshSurfaceApi::cross( const double* a, const double* b, double* c ) const
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];   
}

double FixMeshSurfaceApi::mag( const double* a ) const
{
   
   double res = 0;
   for( int i = 0; i < 3; ++i )
      res += a[i]*a[i];
   return sqrt( res );
   
}
