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

FixStyle(efield/force,FixEfieldForce)

#else

#ifndef LMP_FIX_EFIELD_FORCE_ABSTRACT_H
#define LMP_FIX_EFIELD_FORCE_ABSTRACT_H

#include "fix.h"

#include <cstdlib>
#include "efield_model.h"


namespace LAMMPS_NS {

  class FixEfieldForce : public Fix {
  
  public:
    FixEfieldForce(class LAMMPS *, int, char **);
    ~FixEfieldForce();
    
    virtual int setmask();
    virtual void init();
    virtual void post_force(int vflag);
    
    inline EfieldModel* efieldModel() 
    {
        return efieldModel_;
    } 
    
  protected:
    
    class PairGran *pair_gran;
    
    bool no_torque;
    
    // -- efield model --
    class EfieldModel* efieldModel_;
    
    // -- electric dipole --
    class FixPropertyParticle* fix_p;
    
    // -- particle charge --
    class FixPropertyParticle* fix_charge;

    inline void cross( const double* a, const double* b, double* c ) const
    {
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];   
    }

    
  };

}

#endif
#endif
