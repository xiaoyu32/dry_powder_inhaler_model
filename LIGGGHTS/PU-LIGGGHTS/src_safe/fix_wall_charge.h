
#ifdef FIX_CLASS
FixStyle(wall/gran/efield,FixWallCharge)
#else

#ifndef LMP_FIX_WALL_CHARGE_H
#define LMP_FIX_WALL_CHARGE_H

//#include "contact_interface.h"

#include "fix.h"
#include "fix_charge_gran.h"
#include "fix_charge_gran_conduction.h"
#include "fix_wall_gran.h"
#include <string>
#include <vector>

namespace LCM = LIGGGHTS::ContactModels;

namespace LAMMPS_NS {

class FixWallCharge : public FixWallGran{

public :
  FixWallCharge(class LAMMPS *, int, char **);
  ~FixWallCharge();
 	
  /* INHERITED FROM Fix */
  void post_create();
  void init();
  int setmask();
  void post_force(int vflag);
  
protected :	  

  virtual void compute_force(LCM::CollisionData & cdata, double *vwall);
   
     
  class FixPropertyGlobal* fix_work_wall;
  class FixPropertyGlobal* fix_resistivity_wall;
  class FixEfieldGran * fix_efieldgran;
  class FixPropertyParticle* fix_p;
  
  double transfer_acceleration;
  
  double* work_wall;
  double* resistivity_wall;
  int area_correction_flag;
  int charge_transfer_flag;
  
  double abs( double abs );
  
};

}

#endif
#endif
