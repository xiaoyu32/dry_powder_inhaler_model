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

FixStyle(efield/polarization/permanent,FixPolarizationPermanent)

#else

#ifndef LMP_FIX_POLARIZATION_PERMANENT_H
#define LMP_FIX_POLARIZATION_PERMANENT_H

#include "fix.h"
#include "fix_polarization.h"


namespace LAMMPS_NS {
  
  class FixPolarizationPermanent : public FixPolarization {
    
  friend class FixWallCharge;
  friend class FixEfieldGranCond;  
  friend class FixCfdCouplingElectric;
    
  public:
    FixPolarizationPermanent(class LAMMPS *, int, char **);
    ~FixPolarizationPermanent();
    
    virtual void post_create();
    virtual void pre_delete(bool unfixflag){ UNUSED(unfixflag); };
    
    void initial_integrate(int vflag);
    
    void pre_force(int vflag);
    void final_integrate();
    void end_of_step();
    
    virtual int setmask();
    virtual void init();
    
        
  protected:
    
    double permanent_p[3];
    
  };  

}

#endif
#endif
