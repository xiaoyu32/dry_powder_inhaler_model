/* ----------------------------------------------------------------------
   Jari Kolehmainen and Ali Ozel.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "modify.h"
#include "pointers.h"
#include "math.h"
#include "vector_liggghts.h"
#include "fix_cfd_coupling.h"
//#include "fix_multisphere.h"
#include "cfd_datacoupling_mpi_pun.h"
#include "library_cfd_coupling.h"
#include "properties.h"
#include "modify.h"

#include "particleCommunicator.h"
#include "pushDataMPIpun.h"
#include "pullDataMPIpun.h"



using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CfdDatacouplingMPIpun::CfdDatacouplingMPIpun(LAMMPS *lmp,int iarg, int narg, char **arg,FixCfdCoupling* fc) :
  CfdDatacoupling(lmp, iarg, narg, arg,fc)
{
  
  this->fc_ = fc;
  
  //Modify* modify = this->fc_->modify;
  cpu_mapping = static_cast<FixEulerianCPUs*>(modify->find_fix_id("fcpus"));
    
  /*if(!cpu_mapping)   //FIXME: this does not work, do something!!
  {
  
     printf( "HERE WE ARE!!!!!!!! \n" );
  
    char **fixarg = new char*[3];
    fixarg[0]= (char *) "fcpus";
    fixarg[1]= (char *) "all";
    fixarg[2]= (char *) "cpus";
    fc->modify->add_fix(3,const_cast<char**>(fixarg));
    delete[] fixarg;
    cpu_mapping = static_cast<FixEulerianCPUs*>(fc->modify->find_fix_id("fcpus"));
    
         printf( "HERE WE ARE-2!!!!!!!! \n" );

  }*/   
  
  if( !cpu_mapping ) error->message(FLERR,"CPU mapping is NULL!",1);
    
  liggghts_is_active = false;

  if(!atom->tag_enable) error->one(FLERR,"CFD-DEM coupling via MPI requires particles to have tags");

  if(comm->me == 0) error->message(FLERR,"nevery as specified in LIGGGHTS is overriden by calling external program",1);
  
  particleCommunicator_ = new particleCommunicator(this);
  
  pull_ = new PullDataMPIpun(this); 
  push_ = new PushDataMPIpun(this);
    
  n_particles_CFD = 0;
  
  
}

CfdDatacouplingMPIpun::~CfdDatacouplingMPIpun()
{
   if( pull_ ) delete pull_;
   if( push_ ) delete push_;
   
   pull_ = NULL;
   push_ = NULL;
   
}

void CfdDatacouplingMPIpun::getCurrentCFDSizes()
{
    n_particles_CFD = cpu_mapping->n_cfd;  
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPIpun::exchange()
{
    // does nothing since done by OF
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPIpun::pull(const char *name,const char *type,void *&from,const char *datatype)
{
    
    if( !particleCommunicator_->isUpdated() )
    {
        if(comm->me == 0) printf( "Particle topology updated (during pull call)!\n" );
        particleCommunicator_->do_comm();
    }
    
    CfdDatacoupling::pull(name,type,from,datatype);
    pull_->do_comm( name, type, from, datatype );
    /*
    if(strcmp(datatype,"double") == 0)
        pull_mpi<double>(name,type,from);
    else if(strcmp(datatype,"int") == 0)
        pull_mpi<int>(name,type,from);
    else error->one(FLERR,"Illegal call to CfdDatacouplingMPIpun::pull, valid datatypes are 'int' and double'");*/
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPIpun::push(const char *name,const char *type,void *&to,const char *datatype)
{
    
    // -- should always be
    if( !particleCommunicator_->isUpdated() )
    {
        if(comm->me == 0) printf( "Particle topology updated (during push call)!\n" );
        particleCommunicator_->do_comm();
    }

    CfdDatacoupling::push(name,type,to,datatype);
    push_->do_comm(name,type,to,datatype);
    
    
    /*if(strcmp(datatype,"double") == 0)
        push_mpi<double>(name,type,to);
    else if(strcmp(datatype,"int") == 0)
        push_mpi<int>(name,type,to);
    else error->one(FLERR,"Illegal call to CfdDatacouplingMPIpun::pull, valid datatypes are 'int' and double'");*/
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPIpun::allocate_external(int **&data, int len2,int len1,int initvalue) //FINE
{
  if(len1 < 1 || len2 < 1)
    error->one(FLERR,"Illegal length used in CfdDatacouplingMPI::allocate_external");

  memory->grow(data, len1,len2, "CfdDatacouplingMPI:data");
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
      data[i][j] = initvalue;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPIpun::allocate_external(int    **&data, int len2,char *keyword,int initvalue)
{
  //int len1 = 0;
  //MultisphereParallel *ms_data = properties_->ms_data();

  /*if(strcmp(keyword,"nparticles") == 0) len1 = atom->tag_max();
  else if(strcmp(keyword,"nbodies") == 0)
  {
      if(ms_data)
        len1 = ms_data->tag_max_body();
      else error->one(FLERR,"CFD datacoupling keyword 'nbodies' may only be used with multisphere model in LIGGGHTS");
  }
  else error->one(FLERR,"Illegal keyword used in CfdDatacouplingMPI::allocate_external");
  if(len1 < 1 || len2 < 1)
   len1 = len2 = 1;*/
  CfdDatacouplingMPIpun::getCurrentCFDSizes();
  int len1 = n_particles_CFD;
  
  memory->grow(data, len1,len2, "CfdDatacouplingMPI:data");
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
      data[i][j] = initvalue;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPIpun::allocate_external(double **&data, int len2,int len1,double initvalue)
{
  if(len1 < 1 || len2 < 1)
    error->one(FLERR,"Illegal length used in CfdDatacouplingMPI::allocate_external");
  memory->grow(data, len1,len2, "CfdDatacouplingMPI:data");
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
        data[i][j] = initvalue;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPIpun::allocate_external(double **&data, int len2,char *keyword,double initvalue)
{
  int len1 = 0;
  /*MultisphereParallel *ms_data = properties_->ms_data();

  if(strcmp(keyword,"nparticles") == 0) len1 = atom->tag_max();
  else if(strcmp(keyword,"nbodies") == 0)
  {
      if(ms_data)
        len1 = ms_data->tag_max_body();
      else error->one(FLERR,"CFD datacoupling keyword 'nbodies' may only be used with multisphere model in LIGGGHTS");
  }
  else error->one(FLERR,"Illegal keyword used in CfdDatacouplingMPI::allocate_external");
  if(len1 < 1 || len2 < 1)
    len1 = len2 = 1;*/
  CfdDatacouplingMPIpun::getCurrentCFDSizes();  
  len1 = n_particles_CFD;
  
  memory->grow(data, len1,len2, "CfdDatacouplingMPI:data");
  
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
        data[i][j] = initvalue;
}

void* CfdDatacouplingMPIpun::find_pull_property(const char *name, const char *type, int &len1, int &len2){
    return properties_->find_property_pu(name,type,len1,len2);
}


void* CfdDatacouplingMPIpun::find_push_property(const char *name, const char *type, int &len1, int &len2){
    return properties_->find_property_pu(name,type,len1,len2);
}

/* ---------------------------------------------------------------------- */
