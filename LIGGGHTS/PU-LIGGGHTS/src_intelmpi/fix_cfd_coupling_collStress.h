// Collisional stresses fix class 

#ifdef FIX_CLASS

FixStyle(couple/cfd/collStress,FixCfdCouplingCollStress)

#else

#ifndef LMP_FIX_CFD_COUPLING_COLLSTRESS_H
#define LMP_FIX_CFD_COUPLING_COLLSTRESS_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingCollStress : public Fix  {

  public:
  
  FixCfdCouplingCollStress(class LAMMPS *, int, char **);
  ~FixCfdCouplingCollStress();

  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  virtual void pre_force(int);
  virtual void post_force(int);

  protected:
   
  class FixCfdCoupling* fix_coupling;
  class FixPropertyParticle* fix_collStress;
  
  private:

};

}

#endif
#endif
