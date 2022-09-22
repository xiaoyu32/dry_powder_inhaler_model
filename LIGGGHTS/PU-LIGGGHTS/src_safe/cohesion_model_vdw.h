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
   Yile Gu, Ali Ozel (Sundaresan Group, Princeton)

------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_VDW,vdw,6)
#else
#ifndef COHESION_MODEL_VDW_H_
#define COHESION_MODEL_VDW_H_
#include "contact_models.h"
#include "math.h"

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_VDW> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp), cohEnergyDens(NULL)
    {
      
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("cohEnergyDens", &MODEL_PARAMS::createCohesionEnergyDensity);
      registry.connect("cohEnergyDens", cohEnergyDens,"cohesion_model vdw");
      
      registry.registerProperty("sMin", &MODEL_PARAMS::createSMin);
      registry.connect("sMin", sMin, "cohesion_model vdw");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model vdw");
    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      //r is the distance between the sphere's centers
      //const double r = cdata.r;
      int *type = atom->type;
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const int itype = type[cdata.i];
      const int jtype = type[cdata.j];
      //const double smin = 0.0000400 * (ri + rj); //a minimum cutoff smin of 4E-5 is used to avoit singularity.
      const double smin = sMin[itype][jtype]; //XL: a minimum cutoff smin is used to avoid singularity.
      //const double s = r - (ri + rj);
      //only particle-particle interactions are considered here


      const double Fn_coh = -cohEnergyDens[cdata.itype][cdata.jtype]/3.0
	  * (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin))
	  * ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)   *    ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1) ;

      // Question for normal force: cohesion force should not impact for tangential force, therefore we remove the following statement,
      // cdata.Fn += Fn_coh;

      // apply normal force

        const double fx = Fn_coh * cdata.en[0];
        const double fy = Fn_coh * cdata.en[1];
        const double fz = Fn_coh * cdata.en[2];

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;

        if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    void noCollision(ContactData& cdata, ForceData & i_forces, ForceData & j_forces)
    {
      //double **v = atom->v;
      int *type = atom->type;
      //r is the distance between the sphere's centers
      const double r = sqrt(cdata.rsq);
      const double rinv = 1/r;
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const int itype = type[cdata.i];
      const int jtype = type[cdata.j];
      
      //const double smin = 0.0000400 * (ri + rj); //a minimum cutoff smin of 4E-5 is used to avoit singularity.
      const double smin = sMin[itype][jtype]; //XL: a minimum cutoff smin is used to avoid singularity.
      const double s = r - (ri + rj); // separating distance between the surfaces of the two interacting particles
      //const double smax = (ri + rj) / 4.0; // To speed up the simulation, a maxmimum cutoff separation equal to d/4 is introduced. Beyond smax, the van der Waals cohesive force is not considered.
      const double smax = ri*rj/(ri + rj); // XL
      double Fn_coh;

      if(cdata.touch) *cdata.touch &= ~TOUCH_COHESION_MODEL;

      //only particle-particle interactions are considered here

      if (s < smax)
      {

        if (s < smin)
        {
                Fn_coh = -cohEnergyDens[itype][jtype]/3.0
                * (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin))
                * ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)   *    ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1) ;

        }
        else
        {
                     Fn_coh = -cohEnergyDens[itype][jtype]/3.0
                * (2.0*ri*rj*(ri+rj+s))/(s*s*(2*ri + 2*rj + s)*(2*ri + 2*rj + s))
                * ((s * (2*ri + 2*rj + s)) / ((ri + rj + s)*(ri + rj + s) - (ri-rj)*(ri-rj)) - 1)*((s * (2*ri + 2*rj + s)) / ((ri + rj + s)*(ri + rj + s) - (ri-rj)*(ri-rj)) - 1) ;

        }

        const double fx = Fn_coh * cdata.delta[0] * rinv;
        const double fy = Fn_coh * cdata.delta[1] * rinv;
        const double fz = Fn_coh * cdata.delta[2] * rinv;

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;

        cdata.has_force_update = true;
        if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
    }
  }

  private:
    double ** cohEnergyDens;
    double ** sMin;
  };
}
}

#endif // COHESION_MODEL_VDW_H_
#endif

