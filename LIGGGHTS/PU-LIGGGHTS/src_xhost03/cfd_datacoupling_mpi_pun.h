/* ----------------------------------------------------------------------
  Jari Kolehmainen and Ali Ozel.
------------------------------------------------------------------------- */

#ifdef CFD_DATACOUPLING_CLASS

   CfdDataCouplingStyle(mpipun,CfdDatacouplingMPIpun)

#else

#ifndef LMP_CFD_DATACOUPLING_MPIPUN_H
#define LMP_CFD_DATACOUPLING_MPIPUN_H

#include "pushDataMPIpun.h"
#include "pullDataMPIpun.h"
#include "cfd_datacoupling.h"
#include "particleCommunicator.h"

//#include "multisphere_parallel.h"
#include "error.h"
#include "mpi.h"
#include "fix_eulerian_CPUs.h"
#include "pointers.h"
#include "properties.h"

namespace LAMMPS_NS {

class CfdDatacouplingMPIpun : public CfdDatacoupling {

 friend class PushDataMPIpun;
 friend class PullDataMPIpun;
 friend class particleCommunicator;

 public:
  CfdDatacouplingMPIpun(class LAMMPS *, int,int, char **,class FixCfdCoupling*);
  ~CfdDatacouplingMPIpun();

  void exchange();

  virtual void pull(const char *name, const char *type, void *&ptr, const char *datatype);
  virtual void push(const char *name, const char *type, void *&ptr, const char *datatype);

  template <typename T> void pull_mpi(const char *,const char *,void *&);
  template <typename T> void push_mpi(const char *,const char *,void *&);
  
  inline particleCommunicator* particleComm()
  {
     return particleCommunicator_;
  }
  
  virtual bool error_push()
  { return false;}
    
  void allocate_external(int    **&data, int len2,int len1,     int initvalue);
  void allocate_external(int    **&data, int len2,char *keyword,int initvalue);
  void allocate_external(double **&data, int len2,int len1,     double initvalue);
  void allocate_external(double **&data, int len2,char *keyword,double initvalue);
    
  class FixEulerianCPUs* cpu_mapping;
  
  //non-blocking communicator objects
  class PushDataMPIpun* push_;
  class PullDataMPIpun* pull_;  
  class particleCommunicator* particleCommunicator_;
  
  virtual void markUpdate( bool mark_ )
  {
      particleCommunicator_->setUpdated( mark_ );
  }
  
 protected:
  virtual void* find_pull_property(const char *name, const char *type, int &len1, int &len2);
  virtual void* find_push_property(const char *name, const char *type, int &len1, int &len2);

    
 private:
  
  bool allocate_sizes;
  
  void getCurrentCFDSizes();
  
  int n_particles_CFD;
};

}

#endif
#endif
