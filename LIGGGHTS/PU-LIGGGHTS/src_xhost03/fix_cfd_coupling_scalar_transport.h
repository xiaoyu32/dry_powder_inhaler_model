
#ifdef FIX_CLASS

FixStyle(couple/cfd/scalar_transfer,FixCfdCouplingScalarTransfer)

#else

#ifndef LMP_FIX_CFD_COUPLING_SCALARTRANSFER_H
#define LMP_FIX_CFD_COUPLING_SCALARTRANSFER_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingScalarTransfer : public Fix  {

  public:
  
  FixCfdCouplingScalarTransfer(class LAMMPS *, int, char **);
  ~FixCfdCouplingScalarTransfer();

  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  virtual void pre_force(int);
  virtual void final_integrate();


  class FixCfdCoupling* fix_coupling;
  class FixPropertyParticle** fix_passiveScalar;
  class FixPropertyParticle** fix_passiveScalarFlux;


  protected:
  
  int nPassiveScalars;
  
  double initC;
  bool initValueFlag;
  
  int gran_flag;
  
  class PairGran *pair_gran;
  
};

}

#endif
#endif
