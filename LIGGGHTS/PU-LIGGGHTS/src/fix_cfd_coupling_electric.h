
#ifdef FIX_CLASS

FixStyle(couple/cfd/electric,FixCfdCouplingElectric)

#else

#ifndef LMP_FIX_CFD_COUPLING_ELECTRIC_H
#define LMP_FIX_CFD_COUPLING_ELECTRIC_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingElectric : public Fix  {

  public:
  
  FixCfdCouplingElectric(class LAMMPS *, int, char **);
  ~FixCfdCouplingElectric();

  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  virtual void pre_force(int);
  virtual void post_force(int);
  double compute_vector(int n);

  //ambient electric field coupling
  int electricfield_flag;
  
  int electricfield_force_flag;
  
  //apply correction (used only with the pairwise counting)
  int electicfield_correction;

  class FixCfdCoupling* fix_coupling;
  class FixPropertyParticle* fix_charge;  
  class FixPropertyParticle* fix_electric_field;
  
  
  bool electricfieldGradientFlag;
  class FixPropertyParticle* fix_electric_field_gradient; 
  
  int nullFlag[3];
  
  double constant_ef[3];

  protected:
  
  double ef_total[3];
  
  class FixPropertyParticle* fix_p;
  
  class PairGran *pair_gran;
  
  double permittivity; 
  
  //granular model
  int gran_flag;
  
  private:
  
  double sign( double*, int ) const;
  double sign( double ) const;
  double abs( double ) const;
  inline double dot( double*, double* ) const;
  void cross_product( const double*, const double*, double* ) const;

};

}

#endif
#endif
