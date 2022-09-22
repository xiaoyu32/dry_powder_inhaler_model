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
#include "fix_polarization_dielectric.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPolarizationDielectric::FixPolarizationDielectric(LAMMPS *lmp, int narg, char **arg) : 
FixPolarization(lmp, narg, arg),
permittivity( 0 ),
susceptibility( 0 )
{}

FixPolarizationDielectric::~FixPolarizationDielectric(){}

void FixPolarizationDielectric::post_create()
{
    FixPolarization::post_create();
}

void FixPolarizationDielectric::init()
{
    FixPolarization::init();
    
    // -- connect to properties --
    susceptibility = (static_cast<FixPropertyGlobal*>(modify->find_fix_property("susceptibility","property/global","scalar",0,0,style)))->compute_scalar();
    permittivity = (static_cast<FixPropertyGlobal*>(modify->find_fix_property("permittivity","property/global","scalar",0,0,style)))->compute_scalar();
    
}

int FixPolarizationDielectric::setmask()
{
    return FixPolarization::setmask();
}


void FixPolarizationDielectric::initial_integrate(int vflag)
{    
    int nlocal = atom->nlocal;
    int * mask = atom->mask;


    for( int i = 0; i < nlocal; ++i )
    {
        
	if( !( mask[i] && groupbit ) ) continue; 
        
        for( int j = 0; j < 3; ++j )
	{
	    fix_p->array_atom[i][j] = 0.0;
	    fix_E->array_atom[i][j] = 0.0;
	    fix_Echarge->array_atom[i][j] = 0.0;
	}
    }
    
}

void FixPolarizationDielectric::pre_force(int vflag)
{
    
    //fprintf( screen, "Polarization!!!! %e %e \n", susceptibility, permittivity );
            
    //solve electric field from q
    FixPolarization::update_echarge();

    //iterative solution for p
    double tol = 1.0;
    int it = 0;
    int max_it = 3;
    
    int nlocal = atom->nlocal;
    
    //update_polarization( tol );
            
    while( it < max_it  )
    {
        
	/*for( int i = 0; i < nlocal; ++i )
        {
	    fprintf( screen, "Polarization: %e, %e, %e \n",  
	    fix_p->array_atom[i][0],
	    fix_p->array_atom[i][1],
	    fix_p->array_atom[i][2] );	
	    
	    fprintf( screen, "Efield: %e, %e, %e \n",  
	    fix_Echarge->array_atom[i][0],
	    fix_Echarge->array_atom[i][1],
	    fix_Echarge->array_atom[i][2] );
	    
	}*/
	
        ++it;
        	
	//update particle polarization
        update_polarization( tol );
		
	//update the induced electric field
        update_epolarization();
	
    }
    
}

void FixPolarizationDielectric::end_of_step()
{

}

double FixPolarizationDielectric::compute_scalar()
{
    return 0;
}

void FixPolarizationDielectric::update_polarization( double& tol )
{
    
    int i,j,k;
    int nlocal = atom->nlocal;
    
    int * mask = atom->mask;
    
    double efield[3];
        
    for( i = 0; i < nlocal; ++i )
    {
        
	double res = 0;;
	double res2 = 0;
	
	const double ri = atom->radius[i];
	
	//volume of particle
	const double iV = 4.0*3.1415926*(ri*ri*ri)/3.0;
	
	if( !( mask[i] && groupbit ) ) continue; 
	
	for( k = 0; k < 3; ++k )
	{
	    efield[k] = fix_Echarge->array_atom[i][k]  + fix_E->array_atom[i][k];
	}
	
	for( k = 0; k < 3; ++k )
	{
	    fix_p->array_atom[i][k] = susceptibility * ( iV * permittivity * efield[k] );
	}
	
	
	/*
	fprintf( screen, " fix_Echarge = %e %e %e \n", 	fix_Echarge->array_atom[i][0],
	    								fix_Echarge->array_atom[i][1],
									fix_Echarge->array_atom[i][2] );
									
	fprintf( screen, " fix_E = %e %e %e \n", 	fix_E->array_atom[i][0],
	    								fix_E->array_atom[i][1],
									fix_E->array_atom[i][2] );		
			    
	fprintf( screen, " pvector = %e %e %e \n", 	fix_p->array_atom[i][0],
	    								fix_p->array_atom[i][1],
									fix_p->array_atom[i][2] );
									*/
		

    }
    
    fix_p->do_forward_comm();
}


//integrate the permanent polarization orientation
void FixPolarizationDielectric::final_integrate()
{
    FixPolarization::final_integrate();
}











