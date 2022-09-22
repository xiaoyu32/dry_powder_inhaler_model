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
   VISCOUS DUMMY FORCE MODEL (SIMPLE PHYSICS - NO TESTING)
   Contributing authors:
   Stefan Radl (TU Graz)
   Ali Ozel (Princeton University)
------------------------------------------------------------------------- */



#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_LUBRICATION,lubrication,15)
#else
#ifndef COHESION_MODEL_LUBRICATION_H_
#define COHESION_MODEL_LUBRICATION_H_

#include "contact_models.h"
#include "math.h"
#include "global_properties.h"



namespace LIGGGHTS {

namespace ContactModels {

  template<>
  class CohesionModel<COHESION_LUBRICATION> : protected Pointers {

  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup) :
      Pointers(lmp), liquidViscosity(NULL), historyIndex(0.0), cutOff( 0.001 )
    {
        fprintf( screen, "Lubrication model active.\n" );  
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("liquidViscosity", &MODEL_PARAMS::createCoeffMu);
      registry.connect("liquidViscosity", liquidViscosity,"cohesion_model lubrication"); //connects to FLUID_VISCOSITY (scalar)

       registry.registerProperty("cutOff", &MODEL_PARAMS::createCutOff);      
       registry.connect("cutOff", cutOff,"cohesion_model lubrication");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model lubrication");
	
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

    void cross( const double* a, const double* b, double* c ) const
    {
        
	c[0] = a[1]*b[2] - b[1]*a[2];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
	
    }

// **********************************************************

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      
      if( cdata.is_wall ) return;
      
      double   timeStep = update->dt;
      double        **v = atom->v;
      double**	  omega = atom->omega;
      int         *type = atom->type;
      double      *mass = atom->mass;
      double     *rmass = atom->rmass;
      const double    r = sqrt(cdata.rsq);
      const double radi = cdata.radi;
      const double radj = cdata.radj;
	
	
      if(cdata.touch) *cdata.touch &= ~TOUCH_COHESION_MODEL;

      const double rEff = 2.0 * radi*radj / (radi+radj);

      const double       hatH = cutOff;


      //pull out the viscosity for this pair of particles
      const int itype = type[cdata.i];
      const int jtype = type[cdata.j];
      double myLiqViscosity = liquidViscosity[itype][jtype];

      // relative translational velocity
      const double enx = cdata.delta[0];
      const double eny = cdata.delta[1];
      const double enz = cdata.delta[2];
      
      const double vr1 = v[cdata.i][0] - v[cdata.j][0] ;
      const double vr2 = v[cdata.i][1] - v[cdata.j][1] ;
      const double vr3 = v[cdata.i][2] - v[cdata.j][2] ;
      const double vn  = vr1 * enx + vr2 * eny + vr3 * enz;

      // force scale = force divided by -1*{normal velocity}
      // this is required in order to bound force 
      const double ForceScale = 4.71238898   // 3/2*pi=4.71238898
             * myLiqViscosity 
             * rEff / hatH; 

      const double Fn_coh = ForceScale * -1.0 * vn; //compute force from bounded ForceScale

      //Bound force in order to avoid instability!

      // apply normal force
      const double fx = Fn_coh * cdata.delta[0];
      const double fy = Fn_coh * cdata.delta[1];
      const double fz = Fn_coh * cdata.delta[2];

      i_forces.delta_F[0] += fx;
      i_forces.delta_F[1] += fy;
      i_forces.delta_F[2] += fz;

      j_forces.delta_F[0] -= fx;
      j_forces.delta_F[1] -= fy;
      j_forces.delta_F[2] -= fz;
      
      
      // - tangential lubrication force -
      
      double v_omegati[3];
      double v_omegatj[3];
      
      cross( omega[cdata.i], cdata.delta, v_omegati );
      cross( omega[cdata.j], cdata.delta, v_omegatj );
      
      // - tangential vector
      double nt[3];
      
      double nt_rr = 0;
      
      for( int k = 0; k < 3; ++k )
      {
          v_omegati[k] *= -radi;
	  v_omegatj[k] *= radj; 	  
      }
      
      double vt[3];
      
      vt[0] = vr1 - vn * enx + ( v_omegati[0] - v_omegatj[0] );
      vt[1] = vr2 - vn * eny + ( v_omegati[1] - v_omegatj[1] );
      vt[2] = vr3 - vn * enz + ( v_omegati[2] - v_omegatj[2] );
      
      for( int k = 0; k < 3; ++k )
      {
          nt[k] = vt[k];
	  nt_rr += nt[k] * nt[k];
      }    
      
      nt_rr = sqrt( nt_rr );
      
      if( nt_rr < 1e-8 ) return;
      
      for( int k = 0; k < 3; ++k )
          nt[k] /= nt_rr;
      
      const double vt_mag = nt_rr; //0;
      
      //for( int k = 0; k < 3; ++k )
      //    vt_mag += nt[k] * vt[k];
      
      const double logDinv = hatH < 1.0 ? log( 1.0/ hatH ) : 0;
      
      const double Ft_coh = -vt_mag * 3.1415926 * myLiqViscosity * rEff * logDinv;
      
      const double ftx = Ft_coh * nt[0];
      const double fty = Ft_coh * nt[1];
      const double ftz = Ft_coh * nt[2];     
      
      i_forces.delta_F[0] += ftx;
      i_forces.delta_F[1] += fty;
      i_forces.delta_F[2] += ftz;

      j_forces.delta_F[0] -= ftx;
      j_forces.delta_F[1] -= fty;
      j_forces.delta_F[2] -= ftz;     
      
      //torque;
      double tori[3] = {0,0,0};
      double torj[3] = {0,0,0};
      
      // - direction of torque -
      cross( cdata.delta, nt, tori );
      cross( cdata.delta, nt, torj );
      
      for( int k = 0; k < 3; ++k )
      {
         tori[k] *= radi * Ft_coh;
         torj[k] *= radj * Ft_coh; 
      }
 
      i_forces.delta_torque[0] += tori[0];
      i_forces.delta_torque[1] += tori[1];
      i_forces.delta_torque[2] += tori[2];

      j_forces.delta_torque[0] += torj[0];
      j_forces.delta_torque[1] += torj[1];
      j_forces.delta_torque[2] += torj[2];
      
      // - end of tangential lubrication force -
      
      cdata.has_force_update = true;
      if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;

      cdata.has_force_update = true||cdata.has_force_update;

    }

// **********************************************************
    void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces) {

       double   timeStep = update->dt;
       double        **v = atom->v;
       double**	  omega = atom->omega;
       
       int         *type = atom->type;
       double      *mass = atom->mass;
       double     *rmass = atom->rmass;
       const double    r = sqrt(cdata.rsq);
       const double rinv = 1/r;
       const double radi = cdata.radi;
       const double radj = cdata.radj;
       double dist = r - (radi + radj);

       if( dist < 0 ) dist = 0;	

       if(cdata.touch) *cdata.touch &= ~TOUCH_COHESION_MODEL;


       if ( !cdata.is_wall ) 
       {
         const double rEff = 2.0 * radi*radj / (radi+radj);

         double       hatH = dist / rEff;
	 
	 // - bound distance
	 if( hatH < cutOff ) hatH = cutOff;
	 
	 if( hatH > 1.0 ) return;
	 
         //pull out the viscosity for this pair of particles
         const int itype = type[cdata.i];
         const int jtype = type[cdata.j];
         double myLiqViscosity = liquidViscosity[itype][jtype];

         // relative translational velocity
         const double enx = cdata.delta[0] * rinv;
         const double eny = cdata.delta[1] * rinv;
         const double enz = cdata.delta[2] * rinv;
         const double vr1 = v[cdata.i][0] - v[cdata.j][0];
         const double vr2 = v[cdata.i][1] - v[cdata.j][1];
         const double vr3 = v[cdata.i][2] - v[cdata.j][2];
         const double vn  = vr1 * enx + vr2 * eny + vr3 * enz;

         //force scale = force divided by -1*{normal velocity}
         //this is required in order to bound force 
         const double ForceScale = 4.71238898   // 3/2*pi=4.71238898
        	* myLiqViscosity 
        	* rEff / hatH; 
	
	 const double  Fn_coh = ForceScale * -1.0 * vn ; //compute force from bounded ForceScale

         //Bound force in order to avoid instability!
	
         // apply normal force
         const double fx = Fn_coh * cdata.delta[0] * rinv;
         const double fy = Fn_coh * cdata.delta[1] * rinv;
         const double fz = Fn_coh * cdata.delta[2] * rinv;

         i_forces.delta_F[0] += fx;
         i_forces.delta_F[1] += fy;
         i_forces.delta_F[2] += fz;

         j_forces.delta_F[0] -= fx;
         j_forces.delta_F[1] -= fy;
         j_forces.delta_F[2] -= fz;

	 // - tangential lubrication force -

	 double v_omegati[3];
	 double v_omegatj[3];

	 cross( omega[cdata.i], cdata.delta, v_omegati );
	 cross( omega[cdata.j], cdata.delta, v_omegatj );

	 // - tangential vector
	 double nt[3];

	 double nt_rr = 0;

	 for( int k = 0; k < 3; ++k )
	 {
             v_omegati[k] *= radi * rinv;
	     v_omegatj[k] *= -radj * rinv; 	  
	 }

	 double vt[3];

	 vt[0] = vr1 - vn * enx + ( v_omegati[0] - v_omegatj[0] );
	 vt[1] = vr2 - vn * eny + ( v_omegati[1] - v_omegatj[1] );
	 vt[2] = vr3 - vn * enz + ( v_omegati[2] - v_omegatj[2] );

	 for( int k = 0; k < 3; ++k )
	 {
             nt[k] = vt[k];
	     nt_rr += nt[k] * nt[k];
	 }    

	 nt_rr = sqrt( nt_rr );
	 
	 if( nt_rr < 1e-8 ) return;
	 
	 for( int k = 0; k < 3; ++k )
             nt[k] /= nt_rr;

	 const double& vt_mag = nt_rr;//0;

	 //for( int k = 0; k < 3; ++k )
         //    vt_mag += nt[k] * vt[k];


	 const double logDinv = hatH < 1.0 ? log( 1.0/ hatH ) : 0;
	 
	 //fprintf( screen, "Hejssan dist = %e, reff = %e, hatH = %e \n", dist, rEff, hatH );
	 
	 if( hatH < 1.0 )
	 {
	    	    
	    const double Ft_coh = -vt_mag * 0.5 * 3.1415926 * myLiqViscosity * rEff * logDinv;

	    //fprintf( screen, "n = (%e %e %e) \n", cdata.delta[0]*rinv, cdata.delta[1]*rinv, cdata.delta[2]*rinv );
	    //fprintf( screen, "nt = (%e %e %e) \n", nt[0], nt[1], nt[2] );
	    //fprintf( screen, "vt = (%e %e %e) \n", vt[0], vt[1], vt[2] );
	    //fprintf( screen, "Ft_coh = %e \n", Ft_coh );
	    //fprintf( screen, "%e %e %e %e \n", vt_mag, myLiqViscosity, rEff, logDinv );

	    const double ftx = Ft_coh * nt[0];
	    const double fty = Ft_coh * nt[1];
	    const double ftz = Ft_coh * nt[2];     

	    i_forces.delta_F[0] += ftx;
	    i_forces.delta_F[1] += fty;
	    i_forces.delta_F[2] += ftz;

	    j_forces.delta_F[0] -= ftx;
	    j_forces.delta_F[1] -= fty;
	    j_forces.delta_F[2] -= ftz;     

	    //torque;
	    double tori[3] = {0,0,0};
	    double torj[3] = {0,0,0};

	    // - direction of torque -

	    cross( cdata.delta, nt, tori );
	    cross( cdata.delta, nt, torj );

	    for( int k = 0; k < 3; ++k )
	    {
               tori[k] *= radi * Ft_coh * rinv;
               torj[k] *= radj * Ft_coh * rinv; 
	    }

	    i_forces.delta_torque[0] += tori[0];
	    i_forces.delta_torque[1] += tori[1];
	    i_forces.delta_torque[2] += tori[2];

	    j_forces.delta_torque[0] += torj[0];
	    j_forces.delta_torque[1] += torj[1];
	    j_forces.delta_torque[2] += torj[2];
	
         }
	 
	 // - end of tangential lubrication force -
	  
         cdata.has_force_update = true;
         if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;

       }

       cdata.has_force_update = true||cdata.has_force_update;

    }

    private:

       double ** liquidViscosity;
       double historyIndex;
       double cutOff;

  };
}
}
#endif // COHESION_MODEL_LUBRICATION_H_
#endif


