
#ifdef FIX_CLASS

FixStyle(couple/cfd/api_density,FixCfdCouplingApiDensity)

#else

#ifndef LMP_FIX_CFD_COUPLING_API_TRANSFER_H
#define LMP_FIX_CFD_COUPLING_API_TRANSFER_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingApiDensity : public Fix  {

  public:
  
  FixCfdCouplingApiDensity(class LAMMPS *, int, char **);
  ~FixCfdCouplingApiDensity();

  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  virtual void pre_force(int);
  virtual void final_integrate();


  class FixCfdCoupling* fix_coupling;
  class FixPropertyParticle* fix_api_density;
  class FixPropertyParticle* fix_api_density_collisions_source;
  class FixPropertyParticle* fix_api_density_fluid_source;
  

  protected:
  
  bool fluxCoupled;
    
  double initC;
  bool initValueFlag;
  
  int gran_flag;
  
  class PairGran *pair_gran;
  
};

}

#endif
#endif
