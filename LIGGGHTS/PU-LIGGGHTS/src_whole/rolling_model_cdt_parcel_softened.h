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
#ifdef ROLLING_MODEL
ROLLING_MODEL(ROLLING_CDT_PARCEL_SOFTENED,cdt_parcel_softened,5)
#else
#ifndef ROLLING_MODEL_CDT_PARCEL_SOFTENED_H_
#define ROLLING_MODEL_CDT_PARCEL_SOFTENED_H_
#include "contact_models.h"
#include <algorithm>
#include "math.h"
#include "math_extra_liggghts.h"

namespace LIGGGHTS {
namespace ContactModels
{
  using namespace LAMMPS_NS;

  template<>
  class RollingModel<ROLLING_CDT_PARCEL_SOFTENED> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    RollingModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp), 
	coeffRollFrict(NULL),
        nparcel(NULL), 
	fix_nparcel(NULL)
    {
      
    }

    ~RollingModel()
    {
        if (nparcel) delete [] nparcel;
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("coeffRollFrict", &MODEL_PARAMS::createCoeffRollFrict);
      registry.connect("coeffRollFrict", coeffRollFrict,"rolling_model cdt_parcel_softened");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"rolling model cdt_parcel_softened");
      
      int max_type = atom->get_properties()->max_type();

      fix_nparcel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("nparcel","property/particle","peratomtype",max_type,0,"cdt_parcel_softened"));

      if( !fix_nparcel ) error->all(FLERR,"cdt_parcel_softened requires property/particle fix 'nparcel'!");

      nparcel = new double[max_type];

      for( int i = 0; i < max_type; ++i )
      {
          nparcel[i] = fix_nparcel->compute_vector(i);

          if( nparcel[i] <= 0 ) error->all(FLERR,"Number of particles in a parcel has to be a positive number!");

      }
	
    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces) 
    {
      const double rmu= coeffRollFrict[cdata.itype][cdata.jtype];

      double r_torque[3], wr_roll[3];
      vectorZeroize3D(r_torque);

      // Stiffness scaling
      const double chiStiff = lmp->force->chiStiffnessScaling();

      if(cdata.is_wall){
        const double wr1 = cdata.wr1;
        const double wr2 = cdata.wr2;
        const double wr3 = cdata.wr3;
        const double wrmag = sqrt(wr1*wr1+wr2*wr2+wr3*wr3);
        if (wrmag > 0.)
        {
          const double radius = cdata.radi;
          double kn = cdata.kn;
          const double enx = cdata.en[0];
          const double eny = cdata.en[1];
          const double enz = cdata.en[2];

	  // -- parcel contribution -- 
	  // Fn = kn*deltan, Fn should be multiplied by max(N_i, N_j) (XL)
	  // Since Fn = kn*deltan, we multiplied kn by max(N_i, N_j)
	  kn *= cdata.is_wall? nparcel[cdata.itype-1] : max(nparcel[cdata.itype-1], nparcel[cdata.jtype-1]); 

          r_torque[0] = chiStiff*rmu*kn*(radius-cdata.r)*wr1/wrmag*cdata.cri;
          r_torque[1] = chiStiff*rmu*kn*(radius-cdata.r)*wr2/wrmag*cdata.cri;
          r_torque[2] = chiStiff*rmu*kn*(radius-cdata.r)*wr3/wrmag*cdata.cri;
         
          // remove normal (torsion) part of torque
          double rtorque_dot_delta = r_torque[0]*enx+ r_torque[1]*eny + r_torque[2]*enz;
          double r_torque_n[3];
          r_torque_n[0] = enx * rtorque_dot_delta;
          r_torque_n[1] = eny * rtorque_dot_delta;
          r_torque_n[2] = enz * rtorque_dot_delta;
          vectorSubtract3D(r_torque,r_torque_n,r_torque);
        }
      } else {
        vectorSubtract3D(atom->omega[cdata.i],atom->omega[cdata.j],wr_roll);
        const double wr_rollmag = vectorMag3D(wr_roll);

        if(wr_rollmag > 0.)
        {
          const double radi = cdata.radi;
          const double radj = cdata.radj;
          const double enx = cdata.en[0];
          const double eny = cdata.en[1];
          const double enz = cdata.en[2];
          double kn = cdata.kn;

          // calculate torque
          const double reff= cdata.is_wall ? radi : (radi*radj/(radi+radj));

	  // -- parcel contribution --
	  // Fn = kn*deltan, Fn should be multiplied by max(N_i, N_j) (XL)
	  // Since Fn = kn*deltan, we multiplied kn by max(N_i, N_j)
	  kn *= cdata.is_wall? nparcel[cdata.itype-1] : max(nparcel[cdata.itype-1], nparcel[cdata.jtype-1]); 
          vectorScalarMult3D(wr_roll,chiStiff*rmu*kn*cdata.deltan*reff/wr_rollmag,r_torque);

          // remove normal (torsion) part of torque
          const double rtorque_dot_delta = r_torque[0]*enx + r_torque[1]*eny + r_torque[2]*enz;
          double r_torque_n[3];
          r_torque_n[0] = enx * rtorque_dot_delta;
          r_torque_n[1] = eny * rtorque_dot_delta;
          r_torque_n[2] = enz * rtorque_dot_delta;
          vectorSubtract3D(r_torque,r_torque_n,r_torque);
        }
      }
      
      //r_torque[0] *= cdata.is_wall? nparcel[cdata.itype-1] : max(nparcel[cdata.itype-1], nparcel[cdata.jtype-1]); 
      //r_torque[1] *= cdata.is_wall? nparcel[cdata.itype-1] : max(nparcel[cdata.itype-1], nparcel[cdata.jtype-1]); 
      //r_torque[2] *= cdata.is_wall? nparcel[cdata.itype-1] : max(nparcel[cdata.itype-1], nparcel[cdata.jtype-1]); 
      
      i_forces.delta_torque[0] -= r_torque[0];
      i_forces.delta_torque[1] -= r_torque[1];
      i_forces.delta_torque[2] -= r_torque[2];
      
      // XL: test
      // if (cdata.itype == 1){ std::cout<<"r_torque: torque0: "<<r_torque[0]<<" torque1: "<<r_torque[1]<<" torque2: "<<r_torque[2]<<std::endl; }

      j_forces.delta_torque[0] += r_torque[0];
      j_forces.delta_torque[1] += r_torque[1];
      j_forces.delta_torque[2] += r_torque[2];
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    void noCollision(ContactData&, ForceData&, ForceData&){}

  private:
    double ** coeffRollFrict;

  protected:
      class FixPropertyGlobal* fix_nparcel;
      double * nparcel;
  };
}
}
#endif // ROLLING_MODEL_CDT_PARCEL_SOFTENED_H_
#endif
