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
#include "fix_polarization_permanent.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPolarizationPermanent::FixPolarizationPermanent(LAMMPS *lmp, int narg, char **arg) : FixPolarization(lmp, narg, arg)
{  

   int iarg = 3;
   bool initial_p_flag = false;
   
   permanent_p[0] = 0;
   permanent_p[1] = 0;
   permanent_p[2] = 0;
   
   while( iarg < narg )
   {
       if( strcmp(arg[iarg],"initial_polarization") == 0 )
       {
           ++iarg;
           if( iarg+2 >= narg )
	   {
	       fprintf( screen, "%s \n", arg[iarg-1] ); 
	       error->fix_error(FLERR,this,"not enough arguments for keyword 'initial_polarization'");
	   }

	   permanent_p[0] = atof(arg[iarg]);
	   ++iarg;

	   permanent_p[1] = atof(arg[iarg]);
	   ++iarg;

	   permanent_p[2] = atof(arg[iarg]);

	   initial_p_flag = true;

       }
       ++iarg;
   }
   
   //if( !initial_p_flag ) error->fix_error(FLERR,this,"Specify initial polarization for permanent polarization model");

}

FixPolarizationPermanent::~FixPolarizationPermanent(){}

void FixPolarizationPermanent::post_create()
{

    // -- create polarization vectors with initial value --
    fix_p = static_cast<FixPropertyParticle*>(modify->find_fix_property("polarization","property/particle","vector",3,0,this->style,false));
    if(!fix_p)
    {
      fprintf( screen, "Setting polarization (%e %e %e)...\n", permanent_p[0], permanent_p[1], permanent_p[2] );
    
      char** fixarg = new char*[11];
      fixarg[0]= (char *)"polarization";
      fixarg[1]= (char *)"all";
      fixarg[2]= (char *)"property/particle";
      fixarg[3]= (char *)"polarization";
      fixarg[4]= (char *)"vector";
      fixarg[5]= (char *)"no";
      fixarg[6]= (char *)"yes";
      fixarg[7]= (char *)"no";
      
      fixarg[8] = new char[30];
      fixarg[9] = new char[30];
      fixarg[10] = new char[30];
      
      for( int i = 0; i < 3; ++i )
         sprintf( fixarg[8+i], "%e", permanent_p[i] );
      
      //fixarg[8]= (char *)"0.";
      //fixarg[9]= (char *)"0.";
      //fixarg[10]= (char *)"0.";
      
      fix_p = modify->add_fix_property_atom(11,fixarg,style);
      
      delete[] fixarg[8];
      delete[] fixarg[9];
      delete[] fixarg[10];
      
      delete[] fixarg;
      
    } 

    FixPolarization::post_create();  
}

void FixPolarizationPermanent::init()
{
    FixPolarization::init();
}

int FixPolarizationPermanent::setmask()
{
    return FixPolarization::setmask();
}


void FixPolarizationPermanent::initial_integrate(int vflag)
{    
    FixPolarization::initial_integrate(vflag);
}

void FixPolarizationPermanent::pre_force(int vflag)
{
    
    FixPolarization::pre_force(vflag);
}

void FixPolarizationPermanent::end_of_step()
{
    FixPolarization::end_of_step();
}


//integrate the permanent polarization orientation
void FixPolarizationPermanent::final_integrate()
{
   
   FixPolarization::final_integrate();
   
   int nlocal = atom->nlocal; 
   int*  mask = atom->mask;
   
   double q[4];
   double r;
   double w;
   
   double dt = pair_gran->get_dt();
   
   double** omega = atom->omega;
   double pp,ppo;
   double pi_new[3];
   double* pi_old;
   
   double Rot[3][3];
      
   for( int i = 0; i < nlocal; ++i )
   {
       
       if( !( mask[i] && groupbit ) ) continue;
       
       w = 0;
       pp = 0;
       
       //construct the rotation quaternion
       for( int j = 0; j < 3; ++j )
       {
           pp += fix_p->array_atom[i][j]*fix_p->array_atom[i][j];
           w += omega[i][j]*omega[i][j];   
       }
       
       pi_old = fix_p->array_atom[i];
       
       w = sqrt(w);
       pp = sqrt(pp);
       
       if( w <= 0.0 ) continue;
       
       //rotation quaternion
       q[0] = cos(w/2);
       q[1] = omega[i][0] * sin(w/2)/w;
       q[2] = omega[i][1] * sin(w/2)/w;
       q[3] = omega[i][2] * sin(w/2)/w;
       
       //construct rotation matrix from the quaternion
       Rot[0][0] = q[0]*q[0] + q[1]*q[1] -q[2]*q[2]-q[3]*q[3];
       Rot[0][1] = 2*q[1]*q[2] - 2*q[0]*q[3];
       Rot[0][2] = 2*q[1]*q[3] + 2*q[0]*q[2];
       
       Rot[1][0] = 2*q[1]*q[2] + 2*q[0]*q[3];
       Rot[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
       Rot[1][2] = 2*q[2]*q[3] - 2*q[0]*q[1];
       
       Rot[2][0] = 2*q[1]*q[3] - 2*q[0]*q[2];
       Rot[2][1] = 2*q[2]*q[3] + 2*q[0]*q[1];
       Rot[2][2] = q[0]*q[0] - q[1]*q[1] -q[2]*q[2] + q[3]*q[3];
       
       pi_new[0] = 0;
       pi_new[1] = 0;
       pi_new[2] = 0;
       
       //matrix multiplication
       for( int ii = 0; ii <3; ++ii )
       for( int jj = 0; jj <3; ++jj )
       {
          pi_new[ii] += Rot[ii][jj]*pi_old[jj];
       }
       
       ppo = 0.0;
       
       for( int j = 0; j < 3; ++j )
       {
           ppo += pi_new[j] * pi_new[j];    
       }	
	    
       ppo = sqrt( ppo );	    
       
       if( ppo <= 0.0 ) continue;
       
       //fprintf( screen, "w=%e \n", w );
       //fprintf( screen, "old= %e %e %e \n", fix_permanent_p->array_atom[i][0], fix_permanent_p->array_atom[i][1], fix_permanent_p->array_atom[i][2] );
       //fprintf( screen, "new= %e %e %e \n", pi_new[0], pi_new[1], pi_new[2] );
       //fprintf( screen, "%e %e \n", pp, ppo );
       
       for( int j = 0; j < 3; ++j )
       {
           fix_p->array_atom[i][j] = pi_new[j] * pp/ppo;
       }
       
   }
   
   /*double q[4];
   double r;
   double w;
   
   double dt = pair_gran->get_dt();
   
   double** omega = atom->omega;
   double pp,ppo;
   double pi_new[3];
   double* pi_old;
   
   double Rot[3][3];
   
   int * mask = atom->mask;
   
   for( int i = 0; i < nlocal; ++i )
   {
       
       if( !( mask[i] && groupbit ) ) continue;
       
       w = 0;
       pp = 0;
       
       //relax fixed polarization towards the induced polarization
       for( int j = 0; j < 3; ++j )
       {
       	   fix_permanent_p->array_atom[i][j] += dt/(polar_time)*
	   					( fix_p->array_atom[i][j] - fix_permanent_p->array_atom[i][j] );
       }
       
       
       //construct the rotation quaternion
       for( int j = 0; j < 3; ++j )
       {
           pp += fix_permanent_p->array_atom[i][j]*fix_permanent_p->array_atom[i][j];
           w += omega[i][j]*omega[i][j];   
       }
       
       pi_old = fix_permanent_p->array_atom[i];
       
       w = sqrt(w);
       pp = sqrt(pp);
       
       if( w <= 0.0 ) continue;
       
       //rotation quaternion
       q[0] = cos(w/2);
       q[1] = omega[i][0] * sin(w/2)/w;
       q[2] = omega[i][1] * sin(w/2)/w;
       q[3] = omega[i][2] * sin(w/2)/w;
       
       //construct rotation matrix from the quaternion
       Rot[0][0] = q[0]*q[0] + q[1]*q[1] -q[2]*q[2]-q[3]*q[3];
       Rot[0][1] = 2*q[1]*q[2] - 2*q[0]*q[3];
       Rot[0][2] = 2*q[1]*q[3] + 2*q[0]*q[2];
       
       Rot[1][0] = 2*q[1]*q[2] + 2*q[0]*q[3];
       Rot[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
       Rot[1][2] = 2*q[2]*q[3] - 2*q[0]*q[1];
       
       Rot[2][0] = 2*q[1]*q[3] - 2*q[0]*q[2];
       Rot[2][1] = 2*q[2]*q[3] + 2*q[0]*q[1];
       Rot[2][2] = q[0]*q[0] - q[1]*q[1] -q[2]*q[2] + q[3]*q[3];
       
       pi_new[0] = 0;
       pi_new[1] = 0;
       pi_new[2] = 0;
       
       //matrix multiplication
       for( int ii = 0; ii <3; ++ii )
       for( int jj = 0; jj <3; ++jj )
       {
          pi_new[ii] += Rot[ii][jj]*pi_old[jj];
       }
       
       ppo = 0.0;
       
       for( int j = 0; j < 3; ++j )
       {
           ppo += pi_new[j] * pi_new[j];    
       }	
	    
       ppo = sqrt( ppo );	    
       
       if( ppo <= 0.0 ) continue;
       
       //fprintf( screen, "w=%e \n", w );
       //fprintf( screen, "old= %e %e %e \n", fix_permanent_p->array_atom[i][0], fix_permanent_p->array_atom[i][1], fix_permanent_p->array_atom[i][2] );
       //fprintf( screen, "new= %e %e %e \n", pi_new[0], pi_new[1], pi_new[2] );
       //fprintf( screen, "%e %e \n", pp, ppo );
       
       for( int j = 0; j < 3; ++j )
       {
           fix_permanent_p->array_atom[i][j] = pi_new[j] * pp/ppo;
       }
       
   }*/
    
    
}










