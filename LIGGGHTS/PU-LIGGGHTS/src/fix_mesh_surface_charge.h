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

FixStyle(mesh/surface/charge,FixMeshSurfaceCharge)

#else

#ifndef LMP_FIX_SURFACE_MESH_CHARGE_H
#define LMP_FIX_SURFACE_MESH_CHARGE_H

#include "fix_mesh_surface.h"

namespace LAMMPS_NS
{
  class FixMeshSurfaceCharge : public FixMeshSurface
  {
  
      public:

        FixMeshSurfaceCharge(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMeshSurfaceCharge();

        virtual void post_create();
	
	virtual void compute_scalar_transport( LIGGGHTS::ContactModels::CollisionData& cdata, double * v_wall, int iMesh, int iTri );
	
        virtual void init();
        virtual void setup(int vflag) {}
        virtual int setmask();

        virtual void setup_pre_force(int vflag);

        virtual void pre_force(int vflag);
        virtual void final_integrate();

        virtual double compute_vector(int n);

	// -- compute area of face i --
	double calcArea( int i ) const;
	
	// -- calculate touching point (xMesh) and center to face distance (returned) for particle with center x -- 
	double calcNearestPoint( const double* x, double* xMesh, int iTri ) const;
	
	// -- distance between the face center and x --
	double centerDistance( const double* x, int iTri ) const;
	
	// -- return normal vector pointing from face center to x --
	void normalVector( const double* x, double* n, int iTri ) const;
	
        inline double& surfaceChargeDensity(int i)
        { return (*surfaceChargeDensity_)(i); }

        inline double& surfaceChargeDensityFlux(int i)
        { return (*surfaceChargeDensityFlux_)(i); }
	
	inline const bool& bd_model_flag() const
	{
		return new_bd_model_flag;
	}
	
      protected:
	
	static const double INF = 1e9;
	
	bool chargeTransferActive;
	bool forceActive;
	
	// -- particle charge --
	class FixPropertyParticle* fix_charge;
	class FixPropertyParticle* fix_p;
	
	// -- flux fixes --
	class FixPropertyParticle* fix_efieldFlux;
        class FixPropertyParticle* fix_directionalEfieldFlux;	
	
	// -- previous locations of particles --
	class FixPropertyParticle* fix_prev_loc;
	
	// -- tribocharging model fix --
	class FixEfieldGran* fix_tribo;
	
	TriMesh* triMesh_;
	double initC;
	bool initFlag;
	bool new_bd_model_flag;
		
	// -- wall surface density --
        ScalarContainer<double> *surfaceChargeDensity_;
	
	// -- surface charge density flux --
	ScalarContainer<double> *surfaceChargeDensityFlux_;
        
	// -- (inline) helper functions --	
	inline double distance( const double* a, const double* b ) const
	{
	    double res = 0;
	    for( int i = 0; i < 3; ++i )
	       res += ( a[i] - b[i] ) * ( a[i] - b[i] );
	    return sqrt( res );
	}
	
	
	inline void minus( const double* a, const double* b, double* c ) const
	{
	   for( int i = 0; i < 3; ++i )
	      c[i] = b[i]-a[i];
	}


	inline void cross( const double* a, const double* b, double* c ) const
	{
	    c[0] = a[1]*b[2] - a[2]*b[1];
	    c[1] = a[2]*b[0] - a[0]*b[2];
	    c[2] = a[0]*b[1] - a[1]*b[0];   
	}

	inline double mag( const double* a ) const
	{

	   double res = 0;
	   for( int i = 0; i < 3; ++i )
	      res += a[i]*a[i];
	   return sqrt( res );

	}

	inline double dot( const double* a, const double* b ) const
	{
	    double res = 0;
	    for( int i = 0; i < 3; ++i )
	       res += a[i]*b[i];
	    return res;
	}

	inline void normalize( double* a ) const
	{ 
	    double rr = mag(a);
	    for( int i = 0; i < 3; ++i )
	       a[i] /= rr;
	}

	inline double relu( double a ) const
	{
	    return a>0 ? a : 0;
	}
	

       	// parcel approach
	class FixPropertyGlobal* fix_nparcel;
	double * nparcel;
    
	
  };

} /* namespace LAMMPS_NS */

#endif
#endif
