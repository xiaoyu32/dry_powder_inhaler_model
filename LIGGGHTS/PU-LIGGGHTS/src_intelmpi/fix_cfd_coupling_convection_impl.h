/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2015-     DCS Computing GmbH, Linz
   Copyright 2015-     TU Graz

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

#ifdef FIX_CLASS

FixStyle(couple/cfd/convectiveImpl,FixCfdCouplingConvectiveImpl)

#else

#ifndef LMP_FIX_CFD_COUPLING_CONVECTIVE_IMPL_H
#define LMP_FIX_CFD_COUPLING_CONVECTIVE_IMPL_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingConvectiveImpl : public Fix {

 public:
  FixCfdCouplingConvectiveImpl(class LAMMPS *, int, char **);
  ~FixCfdCouplingConvectiveImpl();
  void post_create();
  void pre_delete(bool unfixflag);

  virtual int setmask();
  virtual void init();
  virtual void post_force(int);

 protected:
  class FixCfdCoupling*     fix_coupling;
  class FixPropertyParticle*    fix_heatFluid;
  class FixPropertyParticle*    fix_heatTransCoeff;
};

}

#endif
#endif
