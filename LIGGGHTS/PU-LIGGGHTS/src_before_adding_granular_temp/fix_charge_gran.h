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

FixStyle(efield/gran,FixEfieldGran)

#else

#ifndef LMP_FIX_EFIELDGRAN_ABSTRACT_H
#define LMP_FIX_EFIELDGRAN_ABSTRACT_H

#include "fix.h"

//Electron charge in coulombs
#define ELECTRON_CHARGE (1.60217657E-19) 

//Change Electron volts to Joules
#define EV2JOULE  (1.60217657E-19)

//define BD_FIELD 48701369.8630137

#define PARTICLE 0
#define WALL 1

#include <cstdlib>
#include "efield_model.h"


namespace LAMMPS_NS {

  class FixEfieldGran : public Fix {
  
  friend class FixInsulatingWall;  
  friend class FixWallCharge;
  friend class FixEfieldGranCond;  
  friend class FixCfdCouplingElectric;
  friend class FixPolarization;
  friend class FixMeshSurfaceChargeConduction;
  friend class FixMeshSurfaceCharge;
  
  public:
    FixEfieldGran(class LAMMPS *, int, char **);
    ~FixEfieldGran();
    
    virtual void post_create();
    virtual void pre_delete(bool unfixflag){ UNUSED(unfixflag); };

    void initial_integrate(int vflag);
    void end_of_step();
    
    virtual double compute_scalar();
    virtual int setmask();
    virtual void init();

    // per default these three methods throw errors.
    virtual void cpl_evaluate(class ComputePairGranLocal *);
    virtual void register_compute_pair_local(class ComputePairGranLocal *);
    virtual void unregister_compute_pair_local(class ComputePairGranLocal *);

    void updatePtrs();
    double copyCharge(int);
    static double work( double r = 1, double a = 0, double b = 0, int type = PARTICLE);
    //int	n_history_extra();
    //bool history_args(char** args);
    double abs( double );
    
    class FixPropertyParticle* fix_ef_coupling;
    
    //neighbor contribution to electricfield 
    class FixPropertyParticle* fix_electricf_r;
    
    double constant_ef[3];
    
    inline const bool& bd_model_flag() const
    {
        return new_bd_model_flag;
    }

    inline const bool& get_relaxation_model_flag() const
    {
	return relaxation_model_flag;
    }
    
    inline EfieldModel* efieldModel() 
    {
        return efieldModel_;
    } 
    
  protected:
  
    // -- efield model --
    EfieldModel* efieldModel_;

  
    class ComputePairGranLocal *cpl;
    class FixPropertyParticle* fix_efieldFlux;
    class FixPropertyParticle* fix_efieldSource;
    class FixPropertyParticle* fix_charge;
    
    class FixScalarTransportEquation *fix_ste;
    
    class FixPropertyParticle* fix_directionalEfieldFlux;
    class FixPropertyParticle* fix_prev_loc;
    
    //Work function parameters
    class FixPropertyGlobal* fix_work_a;
    class FixPropertyGlobal* fix_work_b;
    
    class FixPropertyGlobal* fix_init_charge;
    
    bool random_flag;
    bool opposite_flag;
    
    double* delta_weights;
    double* delta_charges;
    int number_of_delta;
    
    bool delta_flag;
    
    bool multiply_flag;
    bool new_bd_model_flag;
    
    // relaxation model
    bool relaxation_model_flag;
    class FixPropertyGlobal* fix_resistivity;
    double *resistivity;

    double multiply_charge;
    
    bool initChargeOnParticles;
    
    double ** prev_loc;
    double * work_a;
    double * work_b;
    
    double permittivity;
    double delta_charge;
    double bd_field;
    double relaxation_time_scale;
    
    double *efieldFlux;   
    double *efieldSource; 
    double *charge;       
    double C0;          
    double **directionalEfieldFlux;
    bool FHG_init_flag; 
    
    //Global Parameters
    class FixPropertyGlobal* fix_work_model;
    int* work_model;

    class PairGran *pair_gran;
    int dnum,dnum_mine;
    int history_flag;
    int prev_loc_flag;
    
    double random_double();

    int ef_coupling_flag;
    
    //post processing
    class FixPropertyParticle* fix_qflux;
    int flux_flag;
    
  };

}

#endif
#endif
