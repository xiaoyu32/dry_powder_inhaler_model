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
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(mesh/surface/api,FixMeshSurfaceApi)

#else

#ifndef LMP_FIX_SURFACE_MESH_API_H
#define LMP_FIX_SURFACE_MESH_API_H

#include "fix_mesh_surface.h"

namespace LAMMPS_NS
{
  class FixMeshSurfaceApi : public FixMeshSurface
  {
  
      public:

        FixMeshSurfaceApi(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMeshSurfaceApi();

        virtual void post_create();

        virtual void init();
        virtual void setup(int vflag) {}
        virtual int setmask();

        virtual void setup_pre_force(int vflag);

        virtual void pre_force(int vflag);
        virtual void final_integrate();

        virtual double compute_vector(int n);

	// -- compute area of face i --
	double calcArea( int i ) const;
	
      protected:

	TriMesh* triMesh_;
	
	double initC;
	bool initFlag;
	
        inline double& surface_density(int i)
        { return (*surface_density_)(i); }

      private:
	
	// -- wall surface density --
        ScalarContainer<double> *surface_density_;
	
	// -- helper functions --
	void minus( const double* a, const double* b, double* c ) const;
	void cross( const double* a, const double* b, double* c ) const;
	double mag( const double* a ) const;
	
  };

} /* namespace LAMMPS_NS */

#endif
#endif
