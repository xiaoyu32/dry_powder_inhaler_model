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

FixStyle(efield/polarization,FixPolarization)

#else

#ifndef LMP_FIX_POLARIZATION_ABSTRACT_H
#define LMP_FIX_POLARIZATION_ABSTRACT_H

#include "fix.h"
#include "efield_model.h"



namespace LAMMPS_NS {
  
  class FixPolarization : public Fix 
  {
    
  friend class FixWallCharge;
  friend class FixEfieldGranCond;  
  friend class FixCfdCouplingElectric;
  
  friend class EfieldModel;
    
  public:
    FixPolarization(class LAMMPS *, int, char **);
    ~FixPolarization();
    
    virtual void post_create();
    virtual void pre_delete(bool unfixflag){ UNUSED(unfixflag); };
    
    virtual void initial_integrate(int vflag);
    
    virtual void pre_force(int vflag);
    virtual void final_integrate();
    virtual void end_of_step();
    
    virtual double compute_scalar();
    virtual int setmask();
    virtual void init();
        
    void get_electricfield( int, double* );
    
    class FixPropertyParticle* fix_p;
    
  protected:    
    
    EfieldModel* efieldModel;
        
    int ef_coupling_flag;
    
    class FixPropertyParticle* fix_ef_coupling;
    
    class PairGran *pair_gran;
          
    class FixEfieldGran* fix_efield; 
        
    class FixPropertyParticle* fix_E;
    class FixPropertyParticle* fix_Echarge;
    
    void update_echarge();
    void update_epolarization();
    
  };  

}

#endif
#endif
