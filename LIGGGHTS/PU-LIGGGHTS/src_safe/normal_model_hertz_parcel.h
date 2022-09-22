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
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#ifdef NORMAL_MODEL
NORMAL_MODEL(HERTZPARCEL,hertz_parcel,5)
#else
#ifndef NORMAL_MODEL_HERTZ_PARCEL_H_
#define NORMAL_MODEL_HERTZ_PARCEL_H_

#include "fix_property_global.h"
#include "modify.h"
#include "global_properties.h"
#include "math.h"
#include "update.h"


namespace LIGGGHTS {

namespace ContactModels
{
  template<>
  class NormalModel<HERTZPARCEL> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      limitForce(false),
      displayedSettings(false),
      nparcel( NULL ),
      fix_nparcel( NULL )
    {}

    ~NormalModel()
    {
       if( nparcel ) delete[] nparcel;       
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
    }

    void connectToProperties(PropertyRegistry & registry) 
    {
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff,"model hertz");
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff,"model hertz");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model hertz");

      registry.connect("Yeff", Yeff,"model hertz");
      registry.connect("Geff", Geff,"model hertz");
      registry.connect("betaeff", betaeff,"model hertz");

      int max_type = atom->get_properties()->max_type();
      
      fix_nparcel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("nparcel","property/particle","peratomtype",max_type,0,"hertz/parcel"));

      if( !fix_nparcel ) error->all(FLERR,"Normal mode 'hertz_parcel' requires property/particle fix 'nparcel'!");

      nparcel = new double[max_type];

      for( int i = 0; i < max_type; ++i )
      {
	  nparcel[i] = fix_nparcel->compute_vector(i);

	  if( nparcel[i] <= 0 ) error->all(FLERR,"Number of particles in a parcel has to be a positive number!");

      } 
      
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;

      double ri = cdata.radi;
      double rj = cdata.radj;
      double reff = cdata.is_wall ? cdata.radi : (ri*rj/(ri+rj));
      //double meff=cdata.is_wall ? cdata.meff * nparcel[itype-1] : cdata.meff;
      double meff = cdata.meff;

      double sqrtval = sqrt(reff*cdata.deltan);

      double Sn=2.*Yeff[itype][jtype]*sqrtval;
      double St=8.*Geff[itype][jtype]*sqrtval;

      double kn=4./3.*Yeff[itype][jtype]*sqrtval;
      double kt=St;
      const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
      double gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff); 
     
      double gammat=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff);
      
      if (!tangential_damping) gammat = 0.0;

      if(!displayedSettings)
      {
        displayedSettings = true;

        /*
        if(limitForce)
            if(0 == comm->me) fprintf(screen," NormalModel<HERTZ_STIFFNESS>: will limit normal force.\n");
        */
      }
      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*cdata.vn;
      const double Fn_contact = kn*(cdata.radsum-cdata.r);
      const int time = update->ntimestep;
      
      double Fn = Fn_damping + Fn_contact;
           
      // -- parcel contribution --
      // Fn *= cdata.is_wall ? nparcel[itype-1] : nparcel[itype-1] * nparcel[jtype-1];
      Fn *= cdata.is_wall ? nparcel[itype-1] : max(nparcel[itype-1],nparcel[jtype-1]);
      
      //limit force to avoid the artifact of negative repulsion force
      if(limitForce && (Fn<0.0) )
      {
          Fn = 0.0;
      }
      

      cdata.Fn = Fn;
      cdata.kn = kn;
      cdata.kt = kt;
      cdata.gamman = gamman;
      cdata.gammat = gammat;
      
      // apply normal force
      if(cdata.is_wall) {

	const double Fn_ = Fn * cdata.area_ratio;
	
        i_forces.delta_F[0] = Fn_ * cdata.en[0];
        i_forces.delta_F[1] = Fn_ * cdata.en[1];
        i_forces.delta_F[2] = Fn_ * cdata.en[2];

      } else {
        i_forces.delta_F[0] = cdata.Fn * cdata.en[0];
        i_forces.delta_F[1] = cdata.Fn * cdata.en[1];
        i_forces.delta_F[2] = cdata.Fn * cdata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
      }
      
      cdata.has_force_update = true||cdata.has_force_update;
     
      const double maxfraction = 0.25;
      
      // XL:test
      // if (itype == 1) {std::cout<<"Normal Force: Fn1: "<<Fn * cdata.en[0]<<" Fn2: "<<Fn * cdata.en[1]<<" Fn3: "<<Fn * cdata.en[2]<<std::endl; }

      /*
      if(
           cdata.is_wall && 
	   cdata.vn < -1e-16 && 
	   cdata.deltan > maxfraction * reff   
	 ) 
      {
          // -- freeze particles normal velocity --
	  double dt = update->dt;
	  int ii = cdata.i;
	  
	  double cor = exp( betaeff[itype][jtype] * 3.1415926 / ( sqrt( 1.0 -  betaeff[itype][jtype] *  betaeff[itype][jtype] ) ) );
	  
	  if( cor >= 1.0 ) 
	  {
	     cor = 1.0;
	  }

	  if( cor < 0 ) 
	  {
	     cor = 0;
	  }
	  
	  for( int k = 0; k < 3; ++k )
	  { 
	     atom->v[ii][k] = atom->v[ii][k] - (1.0+cor) * cdata.vn * cdata.en[k];
	     //atom->x[ii][k] = atom->x[ii][k] + ( cdata.deltan - maxfraction * reff ) * cdata.en[k];
	  }
	  
	  printf( "Warning: hard sphere wall collision invoked due to exessive wall overlap! (%d) vn =   %e \n", ii, cdata.vn );
	  
      }
      */
      
           
      
    }

    void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** Geff;
    double ** betaeff;
    
    class FixPropertyGlobal* fix_nparcel;
    double * nparcel;
    
    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
  };

}

}
#endif
#endif
