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
COHESION_MODEL(COHESION_VDW_GU,vdw/gu,9)
#else
#ifndef COHESION_MODEL_VDW_GU_H_
#define COHESION_MODEL_VDW_GU_H_
#include "contact_models.h"
#include "math.h"

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_VDW_GU> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp), cohEnergyDens(NULL),
    sMin(NULL), sMins(NULL), sO(NULL)
    {
      if(comm->me==0) printf("COHESION_VDW_GU loaded\n");
      const double chiStiff = lmp->force->chiStiffnessScaling();
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {      
      registry.registerProperty("cohEnergyDens", &MODEL_PARAMS::createCohesionEnergyDensity);
      registry.connect("cohEnergyDens", cohEnergyDens,"cohesion_model vdw/gu");

      registry.registerProperty("sMin", &MODEL_PARAMS::createSMin);
      registry.connect("sMin", sMin, "cohesion_model vdw/gu");
      registry.registerProperty("sMins", &MODEL_PARAMS::createSMins);
      registry.connect("sMins", sMins, "cohesion_model vdw/gu");
      registry.registerProperty("sO", &MODEL_PARAMS::createSO);
      registry.connect("sO", sO, "cohesion_model vdw/gu");  
      
      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model vdw/gu");
	
    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      //r is the distance between the sphere's centers
      const double r = cdata.r;
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const double smin = sMin[cdata.itype][cdata.jtype]; //a minimum cutoff is used to avoit singularity.
      
      // Stifness scaling
      const double chiStiff = lmp->force->chiStiffnessScaling();
      const double realH = cohEnergyDens[cdata.itype][cdata.jtype];
      const double softH = realH/chiStiff;
      
      double Fn_coh;
	if(cdata.is_wall) {
         Fn_coh = -softH * ri /6.0/smin/smin;
	
 	} else {

       	  Fn_coh = -softH/3.0 
		* (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin)) 
		* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)   
		* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1) ;
 	}
     
      // Question for normal force: cohesion force should not impact for tangential force, therefore we remove the following statement,
      // cdata.Fn += Fn_coh;
      
     

      // apply normal force
      if(cdata.is_wall) {
	i_forces.delta_F[0] += Fn_coh * cdata.en[0];
        i_forces.delta_F[1] += Fn_coh * cdata.en[1];
        i_forces.delta_F[2] += Fn_coh * cdata.en[2];
	}
	else {
        const double fx = Fn_coh * cdata.en[0];
        const double fy = Fn_coh * cdata.en[1];
        const double fz = Fn_coh * cdata.en[2];

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;
	}

        if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    void noCollision(ContactData& cdata, ForceData & i_forces, ForceData & j_forces)
    {
      double **v = atom->v;
      int *type = atom->type;
      //r is the distance between the sphere's centers
      const double r = sqrt(cdata.rsq);
      const double rinv = 1/r;
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      double s, smax;	
      if(cdata.is_wall) { 
       	      s = r - ri; // separating distance between the surfaces of the two interacting particles
      	      //smax = ri / 4.0; // To speed up the simulation, a maxmimum cutoff separation equal to d/4 is introduced. Beyond smax, the van der Waals cohesive force is not considered.
	      smax = ri; //XL
      } else{
       	      s = r - (ri + rj); // separating distance between the surfaces of the two interacting particles
       	      //smax = (ri + rj) / 4.0; // To speed up the simulation, a maxmimum cutoff separation equal to d/4 is introduced. Beyond smax, the van der Waals cohesive force is not considered.
	      smax = ri*rj/(ri + rj); //XL
      }	
      double Fn_coh;
                
      //only particle-particle interactions are considered here   
      const int itype = type[cdata.i];
      const int jtype = type[cdata.j];    
      const double smin = sMin[itype][jtype]; //1 nm Yu2012 Chemial Engineering Science
      const double smins = sMins[itype][jtype]; //adjusted smin to take account of softer spring
      const double so = sO[itype][jtype]; //adjusted vdw force equation force to take account of softer spring haha

      // Stifness scaling
      const double chiStiff = lmp->force->chiStiffnessScaling();
      const double realH = cohEnergyDens[itype][jtype];
      const double softH = realH/chiStiff;
      
      if (s < smax)
      {			
	if (s < smins)
	{	
	   if(cdata.is_wall){
             Fn_coh = -softH* ri /6.0/smin/smin;

 	   }else 
	   {
       	     Fn_coh = -softH/3.0 
		     * (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin)) 
		     * ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)   
		     * ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1) ;
	   }              
        }else
	{ 	
	   if(cdata.is_wall) {
             Fn_coh = -realH * ri /6.0/(s-so)/(s-so);
	   } else {
   	     Fn_coh = -realH/3.0 
		     * (2.0*ri*rj*(ri+rj+(s-so)))/((s-so)*(s-so)*(2*ri + 2*rj + (s-so))*(2*ri + 2*rj + (s-so))) 
		     * (((s-so) * (2*ri + 2*rj + (s-so))) / ((ri + rj + (s-so))*(ri + rj + (s-so)) - (ri-rj)*(ri-rj)) - 1)*(((s-so) 
		     * (2*ri + 2*rj + (s-so))) / ((ri + rj + (s-so))*(ri + rj + (s-so)) - (ri-rj)*(ri-rj)) - 1) ;
      	   }
	}
	
    	if(cdata.is_wall) {
	i_forces.delta_F[0] += Fn_coh * cdata.delta[0] * rinv;
        i_forces.delta_F[1] += Fn_coh * cdata.delta[1] * rinv;
        i_forces.delta_F[2] += Fn_coh * cdata.delta[2] * rinv;
	}else 
	{
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
    }  
      

  private:
    double ** cohEnergyDens;
    double ** sMin;
    double ** sMins;
    double ** sO;
  };
}
}


#endif // COHESION_MODEL_VDW_GU_H_
#endif

