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

#ifdef TANGENTIAL_MODEL
TANGENTIAL_MODEL(TANGENTIAL_HISTORY_PARCEL,history_parcel,2)
#else
#ifndef TANGENTIAL_MODEL_HISTORY_PARCEL_H_
#define TANGENTIAL_MODEL_HISTORY_PARCEL_H_
#include "contact_models.h"
#include "math.h"
#include "update.h"
#include "global_properties.h"
#include "atom.h"
#include "fix_property_global.h"
#include "modify.h"

namespace LIGGGHTS {

namespace ContactModels
{
  template<>
  class TangentialModel<TANGENTIAL_HISTORY_PARCEL> : protected Pointers
  {
    double ** coeffFrict;
    int history_offset;

  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    TangentialModel(LAMMPS * lmp, IContactHistorySetup * hsetup) : Pointers(lmp),
      coeffFrict(NULL),
      nparcel( NULL ),
      fix_nparcel( NULL )
    {
      history_offset = hsetup->add_history_value("shearx", "1");
      hsetup->add_history_value("sheary", "1");
      hsetup->add_history_value("shearz", "1");

    }
    
    ~TangentialModel()
    {
    	if ( nparcel ) delete[] nparcel;
    }

    inline void registerSettings(Settings&){}

    inline void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("coeffFrict", &MODEL_PARAMS::createCoeffFrict);
      registry.connect("coeffFrict", coeffFrict,"tangential_model history");
      
      int max_type = atom->get_properties()->max_type();
      
      fix_nparcel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("nparcel","property/particle","peratomtype",max_type,0,"tangential/history/parcel"));

      if( !fix_nparcel ) error->all(FLERR,"history_parcel requires property/particle fix 'nparcel'!");

      nparcel = new double[max_type];

      for( int i = 0; i < max_type; ++i )
      {
	  nparcel[i] = fix_nparcel->compute_vector(i);

	  if( nparcel[i] <= 0 ) error->all(FLERR,"Number of particles in a parcel has to be a positive number!");

      }
      
      
    }

    inline void collision(const CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
        
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      
      // normal forces = Hookian contact + normal velocity damping
      const double enx = cdata.en[0];
      const double eny = cdata.en[1];
      const double enz = cdata.en[2];
            
      // shear history effects
      if(cdata.touch) *cdata.touch |= TOUCH_TANGENTIAL_MODEL;
      double * const shear = &cdata.contact_history[history_offset];
	
      if (cdata.shearupdate && cdata.computeflag) {
        const double dt = update->dt;
        shear[0] += cdata.vtr1 * dt;
        shear[1] += cdata.vtr2 * dt;
        shear[2] += cdata.vtr3 * dt;
	
        // rotate shear displacements

        double rsht = shear[0]*enx + shear[1]*eny + shear[2]*enz;
        shear[0] -= rsht * enx;
        shear[1] -= rsht * eny;
        shear[2] -= rsht * enz;
      }

      const double shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);
      const double kt = cdata.kt;
      const double xmu = coeffFrict[cdata.itype][cdata.jtype];

      // tangential forces = shear + tangential velocity damping
      double Ft1 = -(kt * shear[0]);
      double Ft2 = -(kt * shear[1]);
      double Ft3 = -(kt * shear[2]);

      // rescale frictional displacements and forces if needed
      // const double Ft_shear = kt * shrmag; // sqrt(Ft1 * Ft1 + Ft2 * Ft2 + Ft3 * Ft3);      
      const double Ft_shear = cdata.is_wall ? nparcel[itype-1] * kt * shrmag : max(nparcel[itype-1], nparcel[jtype-1]) * kt * shrmag;
      const double Ft_friction = xmu * fabs(cdata.Fn);

      // energy loss from sliding or damping
      if (Ft_shear > Ft_friction) {
        if (shrmag != 0.0) {
          const double ratio = Ft_friction / Ft_shear;
          Ft1 *= ratio;
          Ft2 *= ratio;
          Ft3 *= ratio;
          shear[0] = -Ft1/kt;
          shear[1] = -Ft2/kt;
          shear[2] = -Ft3/kt;
        }
        else Ft1 = Ft2 = Ft3 = 0.0;
      }
      else
      {
        const double gammat = cdata.gammat;
        Ft1 -= (gammat*cdata.vtr1);
        Ft2 -= (gammat*cdata.vtr2);
        Ft3 -= (gammat*cdata.vtr3);
      }
      
      
      // -- parcel contritbution --
      Ft1 *= cdata.is_wall ? nparcel[itype-1] : max(nparcel[itype-1], nparcel[jtype-1]);
      Ft2 *= cdata.is_wall ? nparcel[itype-1] : max(nparcel[itype-1], nparcel[jtype-1]);
      Ft3 *= cdata.is_wall ? nparcel[itype-1] : max(nparcel[itype-1], nparcel[jtype-1]);

      const double tor1 = eny * Ft3 - enz * Ft2;
      const double tor2 = enz * Ft1 - enx * Ft3;
      const double tor3 = enx * Ft2 - eny * Ft1;
      
      double torque_i[3];
      double torque_j[3];

      torque_i[0] = -cdata.cri * tor1;
      torque_i[1] = -cdata.cri * tor2;
      torque_i[2] = -cdata.cri * tor3;

      if (!cdata.is_wall)
      {   
          torque_j[0] = -cdata.crj * tor1;
          torque_j[1] = -cdata.crj * tor2;
          torque_j[2] = -cdata.crj * tor3;
      }

      // return resulting forces
      if(cdata.is_wall) {
	const double area_ratio = cdata.area_ratio;
        i_forces.delta_F[0] += Ft1 * area_ratio;
        i_forces.delta_F[1] += Ft2 * area_ratio;
        i_forces.delta_F[2] += Ft3 * area_ratio;

        i_forces.delta_torque[0] += torque_i[0] * area_ratio;
        i_forces.delta_torque[1] += torque_i[1] * area_ratio;
        i_forces.delta_torque[2] += torque_i[2] * area_ratio;

      } else {
	i_forces.delta_F[0] += Ft1;
        i_forces.delta_F[1] += Ft2;
        i_forces.delta_F[2] += Ft3;

        j_forces.delta_F[0] -= Ft1;
        j_forces.delta_F[1] -= Ft2;
        j_forces.delta_F[2] -= Ft3;

        i_forces.delta_torque[0] += torque_i[0];
        i_forces.delta_torque[1] += torque_i[1];
        i_forces.delta_torque[2] += torque_i[2];

        j_forces.delta_torque[0] += torque_j[0];
        j_forces.delta_torque[1] += torque_j[1];
        j_forces.delta_torque[2] += torque_j[2];
      }
    }

    inline void noCollision(ContactData & cdata, ForceData&, ForceData&)
    {
      // unset non-touching neighbors
      // TODO even if shearupdate == false?
            
      if(cdata.touch) *cdata.touch &= ~TOUCH_TANGENTIAL_MODEL;
      if(!cdata.contact_history)
          return; //DO NOT access contact_history if not available       
 
      double * const shear = &cdata.contact_history[history_offset];
            
      shear[0] = 0.0;
      shear[1] = 0.0;
      shear[2] = 0.0;
      
    }

    inline void beginPass(CollisionData&, ForceData&, ForceData&){}
    inline void endPass(CollisionData&, ForceData&, ForceData&){}
    
    protected:
      class FixPropertyGlobal* fix_nparcel;
      double * nparcel;
  };
}
}
#endif
#endif
