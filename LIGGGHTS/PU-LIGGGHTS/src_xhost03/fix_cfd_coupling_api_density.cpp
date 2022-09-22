
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include <math.h>
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_api_density.h"
#include "fix_property_atom.h"
#include "atom.h"
#include "compute_pair_gran_local.h"
#include "fix_property_particle.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "mech_param_gran.h"
#include "neigh_list.h"
#include "pair_gran.h"

#define SMALL 1E-16

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingApiDensity::FixCfdCouplingApiDensity(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg)
{
    
    int iarg = 2;
    initValueFlag = false;
    initC = 0;
    gran_flag = 0;
    fluxCoupled = false;
    
    bool hasargs = true;
        
    
    if (strcmp(this->style,"couple/cfd/api_density") != 0) 
    {
	error->fix_error(FLERR,this,"unknown keyword");

    }

    while(iarg < narg && hasargs)
    {
        
	if(strcmp(arg[iarg],"initial_surface_density") == 0)
	{

	    ++iarg;

	    if( iarg >= narg ){
		fprintf( screen, "%s \n", arg[iarg-1] ); 
		error->fix_error(FLERR,this,"not enough arguments for keyword 'initial initial_concentration'");
	    }

	    initC = atof(arg[iarg]);
	    initValueFlag = true;

	}
		
	++iarg;
	
    }
    
    //fprintf( screen, "Number of passive scalars: %d \n", nPassiveScalars );
    
    fix_coupling = NULL;
    fix_api_density = NULL;
    fix_api_density_collisions_source = NULL;
    fix_api_density_fluid_source = NULL;
     
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingApiDensity::~FixCfdCouplingApiDensity()
{
    
    if( fix_api_density ) 			delete fix_api_density;
    if( fix_api_density_collisions_source ) 	delete fix_api_density_collisions_source;
    if( fix_api_density_fluid_source ) 		delete fix_api_density_fluid_source;
    
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingApiDensity::post_create()
{
    
    
    if( !fix_api_density )
    {
	   char* fixarg[9];
	   
	   fixarg[0]= (char *)"carrierConcentration";
           fixarg[1]= (char *)"all";
           fixarg[2]= (char *)"property/particle";
           fixarg[3]= (char *)"carrierConcentration";
           fixarg[4]= (char *)"scalar";
	   fixarg[5]= (char *)"yes"; //restart file
           fixarg[6]= (char *)"yes"; //forward communication with ghost
           fixarg[7]= (char *)"no"; //reverse communication with ghost
	   
	   char arg8[30];
    	   sprintf(arg8,"%e",initC);
           fixarg[8] = arg8;
	   

           fix_api_density = modify->add_fix_property_atom(9,fixarg,style);    
    }
    
    if( !fix_api_density_collisions_source )
    {
	   char* fixarg[9];
	   
	   fixarg[0]= (char *)"carrierCollisionFlux";
           fixarg[1]= (char *)"all";
           fixarg[2]= (char *)"property/particle";
           fixarg[3]= (char *)"carrierCollisionFlux";
           fixarg[4]= (char *)"scalar";
	   fixarg[5]= (char *)"yes"; //restart file
           fixarg[6]= (char *)"no"; //forward communication with ghost
           fixarg[7]= (char *)"no"; //reverse communication with ghost
           fixarg[8] = (char *)"0";
	   
           fix_api_density_collisions_source = modify->add_fix_property_atom(9,fixarg,style);    
    }    

    if( !fix_api_density_fluid_source )
    {
	   char* fixarg[9];
	   
	   fixarg[0]= (char *)"carrierFluidFlux";
           fixarg[1]= (char *)"all";
           fixarg[2]= (char *)"property/particle";
           fixarg[3]= (char *)"carrierFluidFlux";
           fixarg[4]= (char *)"scalar";
	   fixarg[5]= (char *)"yes"; //restart file
           fixarg[6]= (char *)"no"; //forward communication with ghost
           fixarg[7]= (char *)"no"; //reverse communication with ghost
           fixarg[8] = (char *)"0";
	   
           fix_api_density_fluid_source = modify->add_fix_property_atom(9,fixarg,style);    
    }   


}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingApiDensity::pre_delete( bool unfixflag )
{
    if( unfixflag ) modify->delete_fix("fix_api_density");
    if( unfixflag ) modify->delete_fix("fix_api_density_source");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingApiDensity::setmask()
{
    int mask = 0;
    mask |= PRE_FORCE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingApiDensity::init()
{
     
    if(modify->n_fixes_style(style) != 1)
        error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    //check if granular collision model is on, if on employ correction to ambient electric field to
    //avoid double counting of charges
    
    if( force->pair_match("gran", 0) ){ 
        pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
        gran_flag = 1;
    }else{
        gran_flag = 0;
    }

     // find coupling fix
    fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling) error->fix_error(FLERR,this,"Fix couple/cfd/scalar_transfer needs a fix of type couple/cfd");

    if( fix_coupling )
    {
	
	fix_coupling->add_push_property("carrierConcentration", "scalar-atom");
	
	// -- communicating API particles released due to collisions --
	fix_coupling->add_push_property("carrierCollisionFlux", "scalar-atom");
	
	// -- communicating API particles detached/attached due to fluid flow --
	fix_coupling->add_pull_property("carrierFluidFlux", "scalar-atom");
	
    }
    
    
    //printf("Pull Property Electric field added\n");

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingApiDensity::pre_force(int vflag)
{    
     // -- reset the collisional flux --
     if( fluxCoupled )
     {

	 int nlocal = atom->nlocal;
	 double* radius = atom->radius;

	 int *mask = atom->mask;
	 
	 for( int i = 0; i < nlocal; i++ )
	 {
	     if( mask[i] & groupbit ) fix_api_density_collisions_source->vector_atom[i] = 0.0;
	 }
         
	 fluxCoupled = false;
	 
	 int id;
	 MPI_Comm_rank( MPI_COMM_WORLD, &id );
	 
	 if( id == 0 )
	 {
	     fprintf( screen, "Collisional API particle source reseted...\n" );
	 }
	 
     }
     
}

void FixCfdCouplingApiDensity::final_integrate()
{	
    
    if( !fix_api_density )
    {
    	fprintf( screen, "Warning: NULL api density pointer detected!\n" );
        return;
    }
    
    double dt = update->dt;
    int nlocal = atom->nlocal;
    double* radius = atom->radius;
    
    int *mask = atom->mask;
    
    for(int i = 0; i < nlocal; i++)
    {
        if ( mask[i] & groupbit )
	{
	   const double area = 4.0 * M_PI * pow( radius[i], 2 );
           fix_api_density->vector_atom[i] += fix_api_density_fluid_source->vector_atom[i] * dt / area;
	}
    }
    
    // -- mark that collisional flux has been communicated to cfd --
    if( fix_coupling->coupleThis() == 1 )
    {
        fluxCoupled = true;
    }
       
}

/* ----------------------------------------------------------------------
   return components of total force on fix group
------------------------------------------------------------------------- */




















