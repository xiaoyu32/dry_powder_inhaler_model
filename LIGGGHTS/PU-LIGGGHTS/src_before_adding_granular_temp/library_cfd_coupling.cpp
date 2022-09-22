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
   CFD-DEM Coupling Stuff
------------------------------------------------------------------------- */

#include <cstdlib>

#include "mpi.h"
#include "string.h"
#include "library_cfd_coupling.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix_cfd_coupling.h"

//#include "fix_multisphere.h"
#include "cfd_regionmodel.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "cfd_datacoupling.h"
#include "cfd_datacoupling_mpi_pu.h"

#include "fix_mesh.h"
#include "abstract_mesh.h"

using namespace LAMMPS_NS;

#define LMP_GROW_DELTA 11000

/* ---------------------------------------------------------------------- */

bool LocalTagUpdater::wallUpdated_ = true;
bool LocalTagUpdater::flag = false;
int LocalTagUpdater::time_step = 0;
int LocalTagUpdater::localtagCFD_ = 0;


int liggghts_get_maxtag(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return lmp->atom->tag_max();
}

/* ---------------------------------------------------------------------- */

int liggghts_get_maxtag_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies
  
  return 0;
  
  //LAMMPS *lmp = (LAMMPS *) ptr;
  //FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  //if(!fix_ms) return 0;
  //return fix_ms->tag_max_body();
}

/* ---------------------------------------------------------------------- */

int liggghts_get_ntypes_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies
  
  return 0;
  
  //LAMMPS *lmp = (LAMMPS *) ptr;
  //FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  //if(!fix_ms) return 0;
  //return fix_ms->ntypes();
}

/* ---------------------------------------------------------------------- */

double* liggghts_get_vclump_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  
  return 0;
  //FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  //if(!fix_ms) return 0;
  //return fix_ms->vclump();
}

/* ---------------------------------------------------------------------- */

void* locate_coupling_fix(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    int ifix = -1;
    for(int i=0;i<lmp->modify->nfix;i++)
      if(strcmp(lmp->modify->fix[i]->style,"couple/cfd") == 0)
        ifix = i;

    if(ifix ==-1) lmp->error->all(FLERR,"No fix of style 'couple/cfd' found, aborting.");

    return ((void*)lmp->modify->fix[ifix]);
}

/* ---------------------------------------------------------------------- */

void data_liggghts_to_of(char *name,char *type,void *ptr,void *&data,char* datatype)
{
        
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->push(name,type,data,datatype);
}

/* ---------------------------------------------------------------------- */

void data_of_to_liggghts(char *name,char *type,void *ptr,void *data,char* datatype)
{

    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->pull(name,type,data,datatype);
    
}

void wall_data_of_to_liggghts( char* id, double**** fdata, void* lmpPtr )
{
    
    LAMMPS *lmp = (LAMMPS *) lmpPtr;
    
    // -- loop over all the liggghts meshes --
    FixMesh* meshFix = NULL;
    int rank = 0;

    do
    {
       meshFix = static_cast<FixMesh*>( lmp->modify->find_fix_style("mesh/surface", rank) );

       if( meshFix )
       {

	   AbstractMesh* mesh_ = meshFix->mesh();

	   // -- pull data -- 
	   mesh_->pullData( id, (void*&)fdata[rank] );

       }else
       {
	   break;
       }

       ++rank;

    }while( meshFix );
    
}

void wall_data_liggghts_to_of( char* id, double**** tdata, void* lmpPtr  )
{

    LAMMPS *lmp = (LAMMPS *) lmpPtr;
    
    // -- loop over all the liggghts meshes --
    FixMesh* meshFix = NULL;
    int rank = 0;

    do
    {
       meshFix = static_cast<FixMesh*>( lmp->modify->find_fix_style("mesh/surface", rank) );

       if( meshFix )
       {

	   AbstractMesh* mesh_ = meshFix->mesh();

	   // -- pull data -- 
	   mesh_->pushData( id, (void*&)tdata[rank] );

       }else
       {
	   break;
       }

       ++rank;

    }while( meshFix );
   
}

int nMeshes( void* lmpPtr )
{
    LAMMPS *lmp = (LAMMPS *) lmpPtr;
    
    // -- loop over all the liggghts meshes --
    FixMesh* meshFix = NULL;
    int rank = 0;

    do
    {
       meshFix = static_cast<FixMesh*>( lmp->modify->find_fix_style("mesh/surface", rank) );

       if( !meshFix )
       {
	   break;
       }

       ++rank;

    }while( meshFix );
    
    return rank;
    
}

int numberOfElements( void* lmpPtr, int rank )
{
    LAMMPS *lmp = (LAMMPS *) lmpPtr;    
    FixMesh* meshFix = static_cast<FixMesh*>( lmp->modify->find_fix_style("mesh/surface", rank) );
    
    if( !meshFix ) return 0;
    
    AbstractMesh* mesh_ = meshFix->mesh();
    
    if( !mesh_ ) return 0;
    
    //printf( "Number of elements = %d \n", mesh_->ncfd() ); //fixme
    
    return mesh_->ncfd();
    
}

void dataSize( void*& lmpPtr, char* id, int& len, int& nvec )
{
    LAMMPS *lmp = (LAMMPS *) lmpPtr;
    
    // -- loop over all the liggghts meshes --
    FixMesh* meshFix = NULL;
    int rank = 0;
    
    len = 0;
    nvec = 0;
    
    do
    {
       meshFix = static_cast<FixMesh*>( lmp->modify->find_fix_style("mesh/surface", rank) );

       if( meshFix )
       {

	   AbstractMesh* mesh_ = meshFix->mesh();

	   // -- pull data -- 
	   mesh_->dataSize( id, len, nvec );
	   
	   if( len > 0 || nvec > 0 ) return;
	   
       }else
       {
	   break;
       }

       ++rank;

    }while( meshFix );
   
}


/* ---------------------------------------------------------------------- */

void update_rm(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    //FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    locate_coupling_fix(ptr);
    //CfdRegionmodel *rm = fcfd->rm;

    //if(rm) rm->rm_update();
    lmp->error->all(FLERR,"Region model update not implemented aborting.");
}

/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,int len1,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}
/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,char *keyword,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,int len1,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,char* keyword,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

void check_datatransfer(void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->check_datatransfer();
}

/* ---------------------------------------------------------------------- */

int liggghts_get_localtag(void *ptr, bool sanity)
{
   LAMMPS *lmp = (LAMMPS *) ptr;
   int me, nprocs;    
   MPI_Comm_rank(MPI_COMM_WORLD,&me);
   
   //printf( "Processor: %d liggghts_get_localtag  time = %d    tag_time = %d  \n", me, lmp->update->ntimestep, LocalTagUpdater::time_step );
   //MPI_Barrier( MPI_COMM_WORLD );
   
   if( 
       lmp->update->ntimestep == LocalTagUpdater::time_step && 
       LocalTagUpdater::flag 
       && sanity
     )
   { 
     return LocalTagUpdater::localtagCFD_;
   }
    
   if( sanity ) printf( "Error: Sanity check failed! (library_cfd_coupling.cpp, line 223) \n" );//lmp->error->all(FLERR,"Data not yet available!\n");
    
   LocalTagUpdater::time_step = lmp->update->ntimestep;
   LocalTagUpdater::flag = true;

   FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
   CfdDatacouplingMPIpu * cfd_mpi = static_cast<CfdDatacouplingMPIpu*>(fcfd->get_dc());
 

   //MPI_Comm_rank(MPI_COMM_WORLD,&me);
   MPI_Comm_size(MPI_COMM_WORLD,&nprocs); 

   //int n_dem_[nprocs];
   int* n_dem_ = cfd_mpi->cpu_mapping->libraryTagBuffer();
   int* tagCount_ =  cfd_mpi->cpu_mapping->localtagCFD();
   
   // -- initialize to zero --
   for( int i = 0; i < nprocs; ++i )
       n_dem_[i] = 0;
   
   MPI_Allreduce( tagCount_, n_dem_, nprocs, MPI_INT, MPI_SUM, MPI_COMM_WORLD ); //FIXME
   
   LocalTagUpdater::localtagCFD_ = n_dem_[me];   

   //printf( "Finished Reduction! \n", me ); 
   //MPI_Barrier( MPI_COMM_WORLD );
    
   //int n_total = 0;
   //for( int i = 0; i < nprocs; ++i )
   //n_total += n_dem_[i];
   //printf( "n_dem=%d, rank+%d \n", n_dem_[me], me  );      
   return LocalTagUpdater::localtagCFD_;
}

/* ---------------------------------------------------------------------- */

