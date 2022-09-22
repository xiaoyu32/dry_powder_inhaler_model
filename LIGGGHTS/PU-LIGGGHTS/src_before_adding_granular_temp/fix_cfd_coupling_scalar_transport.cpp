
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
#include "fix_cfd_coupling_scalar_transport.h"
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

FixCfdCouplingScalarTransfer::FixCfdCouplingScalarTransfer(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg)
{
    
    int iarg = 2;
    initValueFlag = false;
    initC = 0;
    gran_flag = 0;
    
    bool hasargs = true;
    
    nPassiveScalars = 1;
    
    
   if (strcmp(this->style,"couple/cfd/scalar_transfer") != 0) 
   {
       error->fix_error(FLERR,this,"unknown keyword");

   }
    
    while(iarg < narg && hasargs)
    {
        
	if(strcmp(arg[iarg],"numberOfPassiveScalars") == 0)
	{

	    ++iarg;

	    if( iarg >= narg ){
		fprintf( screen, "%s \n", arg[iarg-1] ); 
		error->fix_error(
				  FLERR,
				  this,
				  "not enough arguments for keyword 'initial numberOfPassiveScalars'"
				);
	    }

	    nPassiveScalars = atoi(arg[iarg]);
	    
	    if( nPassiveScalars < 1 )
		error->fix_error(
				  FLERR,
				  this,
				  "Parameter 'numberOfPassiveScalars' cannot be smaller than 1"
				);	    
	    

	}else if(strcmp(arg[iarg],"initial_concentration") == 0)
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
    fix_passiveScalar = NULL;
    fix_passiveScalarFlux = NULL;
	      
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingScalarTransfer::~FixCfdCouplingScalarTransfer()
{
    
    if( fix_passiveScalar ) free( fix_passiveScalar );
    if( fix_passiveScalarFlux ) free( fix_passiveScalarFlux );
    
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingScalarTransfer::post_create()
{
     
    //if( gran_flag == 0 ) return;

    fprintf( screen, "Number of passive scalars simulated: %d \n", nPassiveScalars );
    
    if( !fix_passiveScalar )
    {
       fix_passiveScalar = (FixPropertyParticle**) malloc( sizeof( FixPropertyParticle* ) * nPassiveScalars );
    }
    
    if( !fix_passiveScalarFlux )
    {
       fix_passiveScalarFlux = (FixPropertyParticle**) malloc( sizeof( FixPropertyParticle* ) * nPassiveScalars );
    }
    
    
    for( int ii = 0; ii < nPassiveScalars; ++ii )
    {
       
       char* strBuff0 = (char*)"passiveScalar";
       char* strBuff1 = (char*)"passiveScalarFlux";
       
       //number of characters needed to present ii
       int nChars( ceil( log10( double(ii+1) ) ) );
        
       
       char* scalarName = (char*) malloc( sizeof( char ) * (13+nChars) );
       char* fluxName = (char*) malloc( sizeof( char ) * (17+nChars) );
       
       for( int jj = 0; jj < 13; ++jj ) scalarName[jj] = strBuff0[jj];

       for( int jj = 0; jj < 17; ++jj ) fluxName[jj] = strBuff1[jj];
       
       sprintf( scalarName+13, "%d", ii );
       sprintf( fluxName+17, "%d", ii );
       
       //fprintf( screen, "Scalar name: %s Flux name: %s !!!!!!!!!!!!!!!! (%d) \n", scalarName, fluxName, ii );
       
       //scalar
       fix_passiveScalar[ii] = static_cast<FixPropertyParticle*>(modify->find_fix_property(scalarName,"property/particle","scalar",0,0,this->style,false));

       if(!fix_passiveScalar[ii])
       {

	   char* fixarg[9];
	   
	   fixarg[0]= scalarName;
           fixarg[1]= (char *)"all";
           fixarg[2]= (char *)"property/particle";
           fixarg[3]= scalarName;
           fixarg[4]= (char *)"scalar";
	   fixarg[5]= (char *)"yes"; //restart file
           fixarg[6]= (char *)"no"; //forward communication with ghost
           fixarg[7]= (char *)"no"; //reverse communication with ghost
	   
	   char arg8[30];
    	   sprintf(arg8,"%e",initC);
           fixarg[8] = arg8;
	   

           fix_passiveScalar[ii] = modify->add_fix_property_atom(9,fixarg,style);

       }
       
       // flux
       fix_passiveScalarFlux[ii] = static_cast<FixPropertyParticle*>(modify->find_fix_property(fluxName,"property/particle","scalar",0,0,this->style,false));

       if(!fix_passiveScalarFlux[ii])
       {

	   char* fixarg[9];

	   fixarg[0]= fluxName;
           fixarg[1]= (char *)"all";
           fixarg[2]= (char *)"property/particle";
           fixarg[3]= fluxName;
           fixarg[4]= (char *)"scalar";           
	   fixarg[5]= (char *)"no";
           fixarg[6]= (char *)"no";
           fixarg[7]= (char *)"no";
	   fixarg[8] = (char*)"0.0";

           fix_passiveScalarFlux[ii] = modify->add_fix_property_atom(9,fixarg,style);

       }
	
       free( scalarName );
       free( fluxName );	
    
    }

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingScalarTransfer::pre_delete( bool unfixflag )
{
    if( unfixflag ) modify->delete_fix("passiveScalar");
    if( unfixflag ) modify->delete_fix("passiveScalarFlux");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingScalarTransfer::setmask()
{
    int mask = 0;
    mask |= PRE_FORCE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingScalarTransfer::init()
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
	
	
	for( int ii = 0; ii < nPassiveScalars; ++ii )
	{
	
	    char* strBuff0 = (char*)"passiveScalar";
	    char* strBuff1 = (char*)"passiveScalarFlux";

	    //number of characters needed to present ii
	    int nChars( ceil( log10( double(ii+1) ) ) );

	    char* scalarName = (char*) malloc( sizeof( char ) * (13+nChars) );
	    char* fluxName = (char*) malloc( sizeof( char ) * (17+nChars) );

	    for( int jj = 0; jj < 13; ++jj ) scalarName[jj] = strBuff0[jj];
	    for( int jj = 0; jj < 17; ++jj ) fluxName[jj] = strBuff1[jj];
	    
	    sprintf( scalarName+13, "%d", ii );
	    sprintf( fluxName+17, "%d", ii );
	    
	    //passive scalar concentration in particles
	    
	    //fprintf( screen, "Scalar Transport Coupling Property %s Added (%d)\n", scalarName, nChars );
	    
	    fix_coupling->add_push_property(scalarName,"scalar-atom");

	    //rate of passive scalar exchanged with fluid

	    fix_coupling->add_pull_property(fluxName,"scalar-atom");
	    
	    free( scalarName );
	    free( fluxName );
	    
	}
	
    }
    
    
    //printf("Pull Property Electric field added\n");

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingScalarTransfer::pre_force(int vflag)
{    

     
}

void FixCfdCouplingScalarTransfer::final_integrate()
{	
    
    if( !fix_passiveScalar ){
    	fprintf( screen, "Warning: NULL passive scalar pointer detected!\n" );
        return;
    }
    double dt = update->dt;
    int nlocal = atom->nlocal;
    
    int *mask = atom->mask;
    
    //update particle passive scalar concentration

    for( int ii = 0; ii < nPassiveScalars; ++ii )
    {

	double* passiveScalar = fix_passiveScalar[ii]->vector_atom;

	for(int i = 0; i < nlocal; i++)
	{
            if ( mask[i] & groupbit )
	    {

		//if( fix_passiveScalarFlux->vector_atom[i] > 0 )	
		//	fprintf( screen, "PassiveScalarFlux: %e \n", fix_passiveScalarFlux->vector_atom[i] );

        	passiveScalar[i] += fix_passiveScalarFlux[ii]->vector_atom[i] * dt;
	    }
	}
	
    }
        
}

/* ----------------------------------------------------------------------
   return components of total force on fix group
------------------------------------------------------------------------- */




















