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

   This contribution was done by:
   Bruno Blais, URPEI, Polytechnique Montreal, Bruno.Blais@polymtl.ca
   Christoph Goniva, DCS Computing GmbH

   it is based on fix_cfd_coupling_force_impilcit and adds an accumulation of the force
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/cfd/force/implicit/accumulated,FixCfdCouplingForceImplicitAccumulated)

#else

#ifndef LMP_FIX_CFD_COUPLING_FORCE_IMPLICIT_ACCUMULATED_H
#define LMP_FIX_CFD_COUPLING_FORCE_IMPLICIT_ACCUMULATED_H

#include "fix_cfd_coupling_force.h"

namespace LAMMPS_NS {

class FixCfdCouplingForceImplicitAccumulated : public FixCfdCouplingForce  {
 public:
  FixCfdCouplingForceImplicitAccumulated(class LAMMPS *, int, char **);
  ~FixCfdCouplingForceImplicitAccumulated();
  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  void post_force(int);
  void end_of_step();

 protected:
  double deltaT_;
  bool   useCN_;
  double CNalpha_;
  class FixPropertyParticle* fix_Ksl_;
  class FixPropertyParticle* fix_uf_;
  class FixPropertyParticle* fix_dragAcc_;
};

}

#endif
#endif
