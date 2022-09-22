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

#include "group.h"
#include "stdlib.h"

#include "update.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_sphere.h"

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
#include "fix_nonspherical.h"
#include "math_vector.h"
#include "domain.h"

#include <fstream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNonSpherical::FixNonSpherical(LAMMPS *lmp, int narg, char **arg) : 
Fix(lmp, narg, arg),
fix_position( NULL ),
fix_normal( NULL ),
fix_velocity( NULL ),
fix_neighbor_id( NULL ),
fix_id( NULL ),
pair_gran( NULL ),
insertFlag( false ),
fix_neighbor_index( NULL )
{
  
  eqdistance = -1;
  kn = -1;
  kt = -1;
  
  gamman = 0;
  gammat = 0;
  
  int iarg = 3;
    
  while( iarg < narg )
  {   
      
      if(strcmp(arg[iarg],"normal_spring_stiffness") == 0)
      {
          ++iarg;
          if( iarg >= narg ){
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'normal_spring_stiffness'");
	  }
	  
	  kn = atof(arg[iarg]); 
	  
	  if( kn <  0 ) error->fix_error(FLERR,this,"Spring stiffness has to be a positive floating point number");
	  
      } else if(strcmp(arg[iarg],"tangential_spring_stiffness") == 0)
      {
          ++iarg;
          if( iarg >= narg ){
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'tangential_spring_stiffness'");
	  }
	  
	  kt = atof(arg[iarg]); 
	  
	  if( kt <  0 ) error->fix_error(FLERR,this,"Spring stiffness has to be a positive floating point number");
      } else if(strcmp(arg[iarg],"eqdistance") == 0)
      {
          ++iarg;
          if( iarg >= narg ){
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'eqdistance'");
	  }
	  
	  eqdistance = atof(arg[iarg]); 
	  
	  if( eqdistance <=  0 ) error->fix_error(FLERR,this,"Particle equilibrium spacing has to be a positive floating point number");
      } else if(strcmp(arg[iarg],"tangential_damper") == 0)
      {
          ++iarg;
          if( iarg >= narg ){
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'tangential_damper'");
	  }
	  
	  gammat = atof(arg[iarg]); 
	  
	  if( gammat <  0 ) error->fix_error(FLERR,this,"Damping coefficient s has to be a positive floating point number");
      } else if(strcmp(arg[iarg],"normal_damper") == 0)
      {
          ++iarg;
          if( iarg >= narg ){
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'normal_damper'");
	  }
	  
	  gamman = atof(arg[iarg]); 
	  
	  if( gammat <  0 ) error->fix_error(FLERR,this,"Damping coefficient s has to be a positive floating point number");
      }    
      
      ++iarg;
  }
  
  if( kn < 0 ) error->fix_error(FLERR,this,"Specify normal spring stiffness using flag 'normal_spring_stiffness'");
  if( kt < 0 ) error->fix_error(FLERR,this,"Specify tangential spring stiffness using flag 'tangential_spring_stiffness'");
  if( eqdistance < 0 ) error->fix_error(FLERR,this,"Specify particle equilibrium spacing using flag 'eqdistance'");
    
}

FixNonSpherical::~FixNonSpherical()
{}

void FixNonSpherical::post_create()
{

  fix_position = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_eqposition","property/particle","vector",3,0,this->style,false));
  if(!fix_position)
  {
  
    char* fixarg[11];
    fixarg[0]= (char *)"nonspherical_eqposition";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"nonspherical_eqposition";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"yes";   // restart
    fixarg[6]= (char *)"yes";   // communicate ghost
    fixarg[7]= (char *)"no";    // communicate rev
    fixarg[8]= (char *)"0.";
    fixarg[9]= (char *)"0.";
    fixarg[10]= (char *)"0.";
    fix_position = modify->add_fix_property_atom(11,fixarg,style);
  } 

  fix_normal = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_normal","property/particle","vector",3,0,this->style,false));
  if(!fix_normal)
  {
  
    char* fixarg[11];
    fixarg[0]= (char *)"nonspherical_normal";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"nonspherical_normal";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"yes";   // restart
    fixarg[6]= (char *)"yes";   // communicate ghost
    fixarg[7]= (char *)"no";    // communicate rev
    fixarg[8]= (char *)"0.";
    fixarg[9]= (char *)"0.";
    fixarg[10]= (char *)"0.";
    fix_normal = modify->add_fix_property_atom(11,fixarg,style);
  } 
  
  fix_velocity = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_velocity","property/particle","vector",3,0,this->style,false));
  if(!fix_velocity)
  {
  
    char* fixarg[11];
    fixarg[0]= (char *)"nonspherical_velocity";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"nonspherical_velocity";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"yes";   // restart
    fixarg[6]= (char *)"yes";   // communicate ghost
    fixarg[7]= (char *)"no";    // communicate rev
    fixarg[8]= (char *)"0.";
    fixarg[9]= (char *)"0.";
    fixarg[10]= (char *)"0.";
    fix_velocity = modify->add_fix_property_atom(11,fixarg,style);
  }   

  fix_neighbor_index = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_index","property/particle","vector",neighborCount,0,this->style,false));
  if(!fix_neighbor_index)
  {
  
    char* fixarg[8+neighborCount];
    fixarg[0]= (char *)"nonspherical_index";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"nonspherical_index";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"yes";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"no";
    
    for( int i = 0; i < neighborCount; ++i )
       fixarg[8+i] = (char*)"-1";
   
    fix_neighbor_index = modify->add_fix_property_atom(8+neighborCount,fixarg,style);
  }
   
  fix_neighbor_id = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_neighbor","property/particle","vector",neighborCount,0,this->style,false));
  if(!fix_neighbor_id)
  {
  
    char* fixarg[8+neighborCount];
    fixarg[0]= (char *)"nonspherical_neighbor";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"nonspherical_neighbor";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"yes";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"no";
    
    for( int i = 0; i < neighborCount; ++i )
       fixarg[8+i] = (char*)"-1";
   
    fix_neighbor_id = modify->add_fix_property_atom(8+neighborCount,fixarg,style);
  }

  fix_id = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_id","property/particle","scalar",1,0,this->style,false));
  if(!fix_id)
  {
  
    char* fixarg[9];    
    fixarg[0]= (char *)"nonspherical_id";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"nonspherical_id";
    fixarg[4]= (char *)"scalar";
    fixarg[5]= (char *)"yes";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"no";
    fixarg[8] = (char*)"0";
   
    fix_id = modify->add_fix_property_atom(9,fixarg,style);
  }  

}

void FixNonSpherical::init()
{
   
   fix_position = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_eqposition","property/particle","vector",3,0,this->style,true));
   fix_normal = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_normal","property/particle","vector",3,0,this->style,true));
   fix_velocity = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_velocity","property/particle","vector",3,0,this->style,true));
   fix_neighbor_index = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_index","property/particle","vector",neighborCount,0,this->style,true));
   
   
   fix_neighbor_id = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_neighbor","property/particle","vector",neighborCount,0,this->style,true));
   fix_id = static_cast<FixPropertyParticle*>(modify->find_fix_property("nonspherical_id","property/particle","scalar",1,0,this->style,true));
   pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
      
}

void FixNonSpherical::setup_pre_exchange()
{
   
   if( !insertFlag )
   {
      initNeighbors();
      insertFlag = true;
   }
}

int FixNonSpherical::setmask()
{
    int mask = 0;
    mask |= PRE_EXCHANGE;
    mask |= PRE_FORCE;
    mask |= POST_FORCE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

void FixNonSpherical::pre_force(int vflag)
{


  int newton_pair = force->newton_pair;
  int inum = pair_gran->list->inum;
  int * ilist = pair_gran->list->ilist;
  int * numneigh = pair_gran->list->numneigh;
  int ** firstneigh = pair_gran->list->firstneigh;
  int *type = atom->type;
  
  double *radius = atom->radius;
  double **x = atom->x;
  double **v = atom->v;
  double **omega = atom->omega;
  
  int nlocal = atom->nlocal;
  int *mask = atom->mask;    
  
  // -- find neighbor indices at the current time step --
  for( int i = 0; i < nlocal; ++i )
  {
     if( mask[i] & groupbit )
     {
         for( int j = 0; j < neighborCount; ++j )
             fix_neighbor_index->array_atom[i][j] = -1.0;
	 
	 fix_position->array_atom[i][0] = 0;
	 fix_position->array_atom[i][1] = 0;
	 fix_position->array_atom[i][2] = 0;
	 
	 fix_normal->array_atom[i][0] = 0;
	 fix_normal->array_atom[i][1] = 0;
	 fix_normal->array_atom[i][2] = 0;

	 fix_velocity->array_atom[i][0] = 0;
	 fix_velocity->array_atom[i][1] = 0;
	 fix_velocity->array_atom[i][2] = 0;
	 
     }   
  }
  
  
  for( int ii = 0; ii < inum; ++ii )
  {
     int i = ilist[ii];
     double xtmp = x[i][0];
     double ytmp = x[i][1];
     double ztmp = x[i][2];
     double radi = radius[i];
     int jnum = numneigh[i];
     int * jlist = firstneigh[i];  
                      	     
     for (int jj = 0; jj < jnum; jj++) 
     {

	int j = jlist[jj];
        j &= NEIGHMASK;
	
	if ( !(mask[i] & groupbit) && !(mask[j] & groupbit) ) continue;
        
	if( !isNeighbor( i, int( fix_id->vector_atom[j] ) ) ) continue;
		
	appendToList( fix_neighbor_index->array_atom[i], j );
	
	// -- add also to ghost particle lists --
	appendToList( fix_neighbor_index->array_atom[j], i );
	 
     }
     
  }
  
  for( int i = 0; i < nlocal; ++i )
  {
     if( mask[i] & groupbit )
     {

         
	 if( listLength( fix_neighbor_index->array_atom[i] ) >= 2 )
	 {
	 	
	    int ii = int( fix_neighbor_index->array_atom[i][0] );
	    int jj = int( fix_neighbor_index->array_atom[i][1] );

            fix_position->array_atom[i][0] = 0.5 * ( x[jj][0]+x[ii][0] );
	    fix_position->array_atom[i][1] = 0.5 * ( x[jj][1]+x[ii][1] );
	    fix_position->array_atom[i][2] = 0.5 * ( x[jj][2]+x[ii][2] );

            fix_velocity->array_atom[i][0] = 0.5 * ( v[jj][0]+v[ii][0] );
	    fix_velocity->array_atom[i][1] = 0.5 * ( v[jj][1]+v[ii][1] );
	    fix_velocity->array_atom[i][2] = 0.5 * ( v[jj][2]+v[ii][2] );
	    
	    fix_normal->array_atom[i][0] = x[jj][0]-x[ii][0];
	    fix_normal->array_atom[i][1] = x[jj][1]-x[ii][1];
	    fix_normal->array_atom[i][2] = x[jj][2]-x[ii][2];
	    
	    normalize( fix_normal->array_atom[i] );
	    
	 }
     }   
  }
     
  // -- communicate position and axis to ghost particles --
  fix_position->do_forward_comm();
  fix_velocity->do_forward_comm();
  fix_normal->do_forward_comm();    
  
}

void FixNonSpherical::post_force(int vflag)
{
  
  int newton_pair = force->newton_pair;
  int inum = pair_gran->list->inum;
  int * ilist = pair_gran->list->ilist;
  int * numneigh = pair_gran->list->numneigh;
  int ** firstneigh = pair_gran->list->firstneigh;
  int *type = atom->type;
  
  double *radius = atom->radius;
  double **x = atom->x;
  double **v = atom->v;
  double **omega = atom->omega;
  
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  
  
  int *mask = atom->mask;    
  //evaluate the external field coming from neighboring particles
  
  //fprintf( screen, "FixNonSpherical::post_force() inum = %d \n", inum );
  
  for( int ii = 0; ii < inum; ++ii )
  {
     int i = ilist[ii];
     double xtmp = x[i][0];
     double ytmp = x[i][1];
     double ztmp = x[i][2];
     double radi = radius[i];
     int jnum = numneigh[i];
     int * jlist = firstneigh[i];  
       
     //fprintf( screen, "FixNonSpherical::post_force() jnum = %d \n", jnum );	     
	     
     for (int jj = 0; jj < jnum; jj++) 
     {
	
	int j = jlist[jj];
        j &= NEIGHMASK;
	
	if ( !(mask[i] & groupbit) && !(mask[j] & groupbit) ) continue;
	
	//fprintf( screen, "Neighbor: i = %d    j = %d     %d \n", int(fix_id->vector_atom[i]), int(fix_id->vector_atom[j]), isNeighbor( i, fix_id->vector_atom[j] ) );
	
	if( !isNeighbor( i, int( fix_id->vector_atom[j] ) ) ) continue;
	
	double dist = 0;
	double vn = 0;
	
	// -- points from i to j --
	double normal[3] = {0,0,0};
	
	for( int kk = 0; kk < 3; ++kk )
	{
	   normal[kk] = x[j][kk] - x[i][kk];  
	   dist += (x[i][kk] -x[j][kk]) * (x[i][kk] -x[j][kk]);
	}
	
	dist = sqrt(dist);
	
	for( int kk = 0; kk < 3; ++kk )
	    normal[kk] /= dist;
	
	for( int kk = 0; kk < 3; ++kk )
	    vn += normal[kk] * ( v[j][kk] - v[i][kk] );
	
	// -- difference from equilibrium spacing --
	const double delta = dist - eqdistance;
	
	// normal force
	double fn[3] = {0,0,0};
	
	// -- normal force (linear spring and damper) --
	for( int kk = 0; kk < 3; ++kk )
	   fn[kk] = normal[kk] * ( kn * delta + gamman * vn );
		
	for( int kk = 0; kk < 3; ++kk )
	{
	   atom->f[i][kk] += fn[kk];
	}
	
	
	if( j < nlocal )
	{
            for( int kk = 0; kk < 3; ++kk )
	    {
	       atom->f[j][kk] -= fn[kk];
	    }	
	}	

     }
  }
  
  
  // -- tangential forces using terniary interactions --
  for( int i = 0; i < nlocal + nghost; ++i )
  {
      
      double ft[3] = {0,0,0};
      
      double* xeq = fix_position->array_atom[i];
      double* veq = fix_velocity->array_atom[i];
      
      
      double* axis = fix_normal->array_atom[i]; 
       
      //fprintf( screen, "axis = [ %e %e %e ]\n", axis[0], axis[1], axis[2] ); 
       
      if( mag(axis) < 1e-10 ) continue; 
       
      double dummy[3] = {0,0,0}; 
      double normal[3] = {0,0,0};
      
      double del[3] = {0,0,0};
      del[0] = x[i][0] - xeq[0];
      del[1] = x[i][1] - xeq[1];
      del[2] = x[i][2] - xeq[2];
      
      if( mag( del ) < 1e-14 ) continue;
      
      double delv[3] = {0,0,0};
      delv[0] = v[i][0] - veq[0];
      delv[1] = v[i][1] - veq[1];
      delv[2] = v[i][2] - veq[2];  
  
  
      cross( axis, del, dummy );
      normalize( dummy );
      
      // -- compute normal vector between particle position and the line segment between neighboring particles --
      cross( dummy, axis, normal );
      normalize( normal );
      
      //fprintf( screen, "\n i = %d \n", i ); 
      
      //fprintf( screen, "normal = [%e %e %e] \n", normal[0], normal[1], normal[2] );
      //fprintf( screen, "del = [%e %e %e] \n", del[0], del[1], del[2] );
      
      double delta = abs( dot( del, normal ) );
      
      //fprintf( screen, "delta = %e \n", delta );
      
      double vt = dot( delv, axis );
      
      double normalV[3] = {0,0,0};
      
      normalV[0] = delv[0] - vt * axis[0];
      normalV[1] = delv[1] - vt * axis[1];
      normalV[2] = delv[2] - vt * axis[2];

      //fprintf( screen, "normalV = [%e %e %e] \n", normalV[0], normalV[1], normalV[2] );      
      
      for( int kk = 0; kk < 3; ++kk )
          ft[kk] = kt * delta * normal[kk] + gammat * normalV[kk];

      if( i < nlocal )
      {
          for( int kk = 0; kk < 3; ++kk )
	  {
	     atom->f[i][kk] -= ft[kk];
	  }	
      }

      for( int j = 0; j < neighborCount; ++j )
      {
          
	  int index = int( fix_neighbor_index->array_atom[i][j] );
	  
	  if( index < 0 ) break;
	  
	  if( index < nlocal )
	  {
              for( int kk = 0; kk < 3; ++kk )
	      {
		 atom->f[index][kk] += ft[kk]/neighborCount;
	      }	      
	  }
	  
      } 
      
  }
  
}

void FixNonSpherical::final_integrate()
{
}

void FixNonSpherical::initNeighbors()
{
    
    int npart, me;
    
    MPI_Comm_rank( MPI_COMM_WORLD, &me );
    
    // -- starting point --
    double* xpart;
    
    // -- particle direction --
    double* normals;
    
    // -- number of particles in a single non-spherical particle --
    int* mpart;
     
    // -- diameter of spheres --
    double* dpart; 
    
    double* partDensity;
     
    // -- particle type -- 
    int* partType; 
     
    int nlocal_old = atom->nlocal;

    // clear ghost count and any ghost bonus data internal to AtomVec
    // same logic as beginning of Comm::exchange()
    // do it now b/c creating atoms will overwrite ghost atoms
 
    atom->nghost = 0;
    atom->avec->clear_bonus();
     
    // -- read non spherical particles at root --
    if( me == 0 )
    {
    
	std::ifstream file( "nonspherical_particles.dat" );    

	file>>npart;
	
	if( npart > 0 ) 
	{
	   xpart =new double[npart*3];
	   normals = new double[npart*3];
	   mpart = new int[npart];

	   dpart =new double[npart];
	   partDensity =new double[npart];
	   partType = new int[npart];


	   for( int i = 0; i < npart; ++i )
	   {
	       for( int j = 0; j < 3; ++j )
	           file>>xpart[3*i+j];    

	       file>>mpart[i];
	       file>>dpart[i];
	       file>>partDensity[i];
	       file>>partType[i];

	       for( int j = 0; j < 3; ++j )
	           file>>normals[3*i+j];
	   }
	}
    }
    
    // -- send partilce information to other MPI tasks --
    MPI_Bcast( &npart, 1, MPI_INT, 0, MPI_COMM_WORLD );
    if( npart == 0 ) return;
    
    int idcounter = 0;
    
    double currentXpart[3] = {0,0,0};
    double currentNormal[3] = {0,0,0};
    
    double currentDpart = 0;
    double currentDensity = 0;
    
    int currentType = 0;
    int currentMpart = 0;
    
    int nfix = modify->nfix;
    Fix **fix = modify->fix;
    
    // -- insert particles --
    for( int i = 0; i < npart; ++i )
    {
        
	if( me == 0 )
	{
	    currentXpart[0] = xpart[3*i];
	    currentXpart[1] = xpart[3*i+1];
	    currentXpart[2] = xpart[3*i+2];

	    currentNormal[0] = normals[3*i];
	    currentNormal[1] = normals[3*i+1];
	    currentNormal[2] = normals[3*i+2];	    
	    
	    currentType = partType[i];
	    currentDpart = dpart[i];
	    currentDensity = partDensity[i];
	    
	    currentMpart = mpart[i];
	}
	
	// -- bcast only one particle at a time to avoid memory over flow --
	MPI_Bcast( &currentMpart, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &currentType, 1, MPI_INT, 0, MPI_COMM_WORLD );
	    
	MPI_Bcast( &currentDpart, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &currentDensity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	
	MPI_Bcast( &currentXpart, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );    
	MPI_Bcast( &currentNormal, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );     
	    
	for( int j = 0; j <= currentMpart; ++j )
	{ 
	   double xi[3];
	
	   xi[0] = currentXpart[0] + currentDpart * j * currentNormal[0];
	   xi[1] = currentXpart[1] + currentDpart * j * currentNormal[1];
	   xi[2] = currentXpart[2] + currentDpart * j * currentNormal[2];
	
	   if( domain->is_in_subdomain( xi ) )
	   {
	        atom->avec->create_atom(currentType,xi);

                int m = atom->nlocal - 1;
                atom->mask[m] = 1 | groupbit;
		atom->v[m][0] = 0;
                atom->v[m][1] = 0;
		atom->v[m][2] = 0;
		
		atom->omega[m][0] = 0;
		atom->omega[m][1] = 0;
		atom->omega[m][2] = 0;
		
                atom->radius[m] = currentDpart/2.0;
                atom->density[m] = currentDensity;
                atom->rmass[m] = currentDensity * ( 3.1415926 * pow(currentDpart,3) / 6.0 );
		
		//atom->tag[i] = idcounter;
		
		atom->data_fix_compute_variable(atom->nlocal-1,atom->nlocal);
		
                //pre_set_arrays() called above
                for (int k = 0; k < nfix; ++k )
                   if (fix[k]->create_attribute) fix[k]->set_arrays(m);
		   
		// -- set non-spherical particle id, neighbors, and direction --
		fix_id->vector_atom[m] = double( idcounter );   
		
		int neighborCounter = 0;
		
		if( j == 0 ) // -- tail -- (only one neighbor)
		    fix_neighbor_id->array_atom[m][neighborCounter++] = idcounter+1;
		else if( j == currentMpart ) // -- head -- (only one neighbor)	
		    fix_neighbor_id->array_atom[m][neighborCounter++] = idcounter-1;
		else
		{
		    fix_neighbor_id->array_atom[m][neighborCounter++] = idcounter+1;
		    fix_neighbor_id->array_atom[m][neighborCounter++] = idcounter-1;
		} 
		
		//nonspherical_eqposition->array_atom[m][0] = 0;
		//nonspherical_eqposition->array_atom[m][1] = 0;
		//nonspherical_eqposition->array_atom[m][2] = 0;
		
		for( int k = neighborCounter; k < neighborCount; ++k )
		    fix_neighbor_id->array_atom[m][k] = -1;
		  
	       
	   }

	   ++idcounter;
	
	}
	
    }
    
    MPI_Barrier( MPI_COMM_WORLD );
    
    if (atom->tag_enable) atom->tag_extend();
    atom->tag_check();
    
    // if global map exists, reset it
    // invoke map_init() b/c atom count has grown
    if (atom->map_style) 
    {
      atom->map_init();
      atom->map_set();
    }

    if( me == 0 )
    {
       fprintf( screen, "Inserted %d non-spherical particles using %d spheres.\n", npart, idcounter );
              
       delete[] xpart;
       delete[] normals;
       delete[] mpart;
       
       delete[] dpart;
       delete[] partDensity;
       delete[] partType;

    }
    
}

bool FixNonSpherical::isNeighbor( int ilocal, int id ) const
{
    
    for( int ii = 0; ii < neighborCount; ++ii )
    {
        if( int(fix_neighbor_id->array_atom[ilocal][ii]) == id ) return true;
    }
    
    return false;
    
}


