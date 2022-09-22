/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Efield
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(efield/gran/particle,FixEfieldGranCond)

#else

#ifndef LMP_FIX_EFIELDGRAN_CONDUCTION_H
#define LMP_FIX_EFIELDGRAN_CONDUCTION_H

#include "fix_charge_gran.h"
#include <cstdlib>

namespace LAMMPS_NS {

  class FixEfieldGranCond : public FixEfieldGran {
    
  public:
  
    FixEfieldGranCond(class LAMMPS *, int, char **);
    ~FixEfieldGranCond();	
    void pre_delete(bool);

    int setmask();
    void init();
    //void initial_integrate(int vflag);

    void post_force(int);

    void cpl_evaluate(class ComputePairGranLocal *);
    void register_compute_pair_local(ComputePairGranLocal *);
    void unregister_compute_pair_local(ComputePairGranLocal *);
    
    class FixPropertyParticle* fix_p;
     
  private:
    template <int> void post_force_eval(int,int);
    
    bool polarization_flag;
    
    // for efield transfer area correction
    int area_correction_flag;
    double transfer_acceleration;    
    double const* const* deltan_ratio;
    
    
    
  };

}

#endif
#endif

