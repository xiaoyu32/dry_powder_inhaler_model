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

/* ----------------------------------------------------------------------
   Charge transfer between insulating walls and particles
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(mesh/surface/charge/conduction,FixMeshSurfaceChargeConduction)

#else

#ifndef LMP_FIX_SURFACE_MESH_CHARGE_CONDUCTION_H
#define LMP_FIX_SURFACE_MESH_CHARGE_CONDUCTION_H

#include "contact_interface.h"
#include "fix_mesh_surface_charge.h"

namespace LAMMPS_NS
{
  class FixMeshSurfaceChargeConduction : public FixMeshSurfaceCharge
  {
  
      public:

        FixMeshSurfaceChargeConduction(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMeshSurfaceChargeConduction();

	virtual void compute_scalar_transport( LIGGGHTS::ContactModels::CollisionData& cdata, double * v_wall, int iMesh, int iTri );
	
        virtual void post_create();
        virtual void init();
	virtual int setmask();
	
	inline const double& normalElectricField( int i ) const
	{ return (*normalElectricField_)(i); }
	
      protected:
	class FixPropertyGlobal* fix_acc;
	class FixPropertyGlobal* fix_work_wall;
	class FixPropertyGlobal* fix_resistivity_wall;
		
	class FixPropertyGlobal* fix_work_model;
        class FixPropertyGlobal* fix_work_a;
        class FixPropertyGlobal* fix_work_b;	
	
        class FixPropertyParticle* fix_ef_coupling; // -- testing --

	ScalarContainer<double>* normalElectricField_;
	double transfer_acceleration;

  };

} /* namespace LAMMPS_NS */

#endif
#endif
