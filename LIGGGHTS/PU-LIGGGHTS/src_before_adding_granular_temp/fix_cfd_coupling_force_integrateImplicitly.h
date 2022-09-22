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
------------------------------------------------------------------------- */
/* Description
   This fix will couple to the CFD provided forces, and integrate the
   particle position and velocity considering forces in f[i] explicitly,
   and drag forces (provided via coefficient Ksl) IMPLICITLY.
   NO additional fix for integration is necessary!
   This fix always uses a Crank-Nicholson style integration
   for the drag foces (i.e., blending of old and new particle velocity
   based on CNalpha_)
*/

#ifdef FIX_CLASS

FixStyle(couple/cfd/force/integrateImp,FixCfdCouplingForceIntegrateImplicitly)

#else

#ifndef LMP_FIX_CFD_COUPLING_FORCE_INTEGRATEMPLICITLY_H
#define LMP_FIX_CFD_COUPLING_FORCE_INTEGRATEMPLICITLY_H

#include "fix_cfd_coupling_force.h"
#include "fix_nve.h"

namespace LAMMPS_NS {

class FixCfdCouplingForceIntegrateImplicitly : public FixCfdCouplingForce {
 public:
  FixCfdCouplingForceIntegrateImplicitly(class LAMMPS *, int, char **);
  ~FixCfdCouplingForceIntegrateImplicitly();
          void post_create();
          void pre_delete(bool unfixflag);

          int  setmask();
  virtual void init();
  virtual void post_force(int);

  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void initial_integrate_respa(int, int, int);
  virtual void final_integrate_respa(int, int);

          void end_of_step();
  virtual void reset_dt();

 protected:
  double dtv,dtf;
  double *step_respa;
  int extra;

  double CNalpha_;  //Crank-Nicholson blending factor

  double CAddRhoFluid_;   //Added mass coefficient times relative fluid density (C_add*rhoFluid/rhoP)
  double onePlusCAddRhoFluid_;

  class FixPropertyParticle* fix_Ksl_;
  class FixPropertyParticle* fix_uf_;

  inline void implicitVelocityUpdate
            (
                double dtf, double rmass,
                double *v, double *f, double Ksl, double *uf,
                double *frc
            );

};

}

#endif
#endif
