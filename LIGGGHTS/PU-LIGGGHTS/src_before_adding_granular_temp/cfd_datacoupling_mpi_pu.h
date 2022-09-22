/* ----------------------------------------------------------------------
  Jari Kolehmainen and Ali Ozel.
------------------------------------------------------------------------- */

#ifdef CFD_DATACOUPLING_CLASS

   CfdDataCouplingStyle(mpipu,CfdDatacouplingMPIpu)

#else

#ifndef LMP_CFD_DATACOUPLING_MPIPU_H
#define LMP_CFD_DATACOUPLING_MPIPU_H

#include "pushDataMPIpu.h"
#include "pullDataMPIpu.h"
#include "cfd_datacoupling.h"
//#include "multisphere_parallel.h"
#include "error.h"
#include "mpi.h"
#include "fix_eulerian_CPUs.h"
#include "pointers.h"
#include "properties.h"

namespace LAMMPS_NS {

class CfdDatacouplingMPIpu : public CfdDatacoupling {

 friend class PushDataMPIpu;
 friend class PullDataMPIpu;

 public:
  CfdDatacouplingMPIpu(class LAMMPS *, int,int, char **,class FixCfdCoupling*);
  ~CfdDatacouplingMPIpu();

  void exchange();

  virtual void pull(const char *name, const char *type, void *&ptr, const char *datatype);
  virtual void push(const char *name, const char *type, void *&ptr, const char *datatype);

  template <typename T> void pull_mpi(const char *,const char *,void *&);
  template <typename T> void push_mpi(const char *,const char *,void *&);

  virtual bool error_push()
  { return false;}
    
  void allocate_external(int    **&data, int len2,int len1,     int initvalue);
  void allocate_external(int    **&data, int len2,char *keyword,int initvalue);
  void allocate_external(double **&data, int len2,int len1,     double initvalue);
  void allocate_external(double **&data, int len2,char *keyword,double initvalue);
    
  class FixEulerianCPUs* cpu_mapping;
  
  class PushDataMPIpu* push_;
  class PullDataMPIpu* pull_;  

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
