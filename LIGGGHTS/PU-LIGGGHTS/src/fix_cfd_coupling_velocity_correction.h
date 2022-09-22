
// Velocity correction fix class for periodic flow configuration

#ifdef FIX_CLASS

FixStyle(couple/cfd/velocityCorrection,FixCfdCouplingVelocityCorrection)

#else

#ifndef LMP_FIX_CFD_COUPLING_VELOCITY_CORRECTION_H
#define LMP_FIX_CFD_COUPLING_VELOCITY_CORRECTION_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingVelocityCorrection : public Fix  {

  public:
  
  FixCfdCouplingVelocityCorrection(class LAMMPS *, int, char **);
  ~FixCfdCouplingVelocityCorrection();

  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  virtual void pre_force(int);
  virtual void post_force(int);
  double compute_vector(int n);

  protected:
  
  double velCorr_total[3];
  
  class FixCfdCoupling* fix_coupling;
  class FixPropertyParticle* fix_velocity_correction;
  
  private:

};

}

#endif
#endif
