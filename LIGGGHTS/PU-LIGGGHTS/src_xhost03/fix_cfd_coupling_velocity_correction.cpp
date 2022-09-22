// Velocity correction fix class for periodic flow configuration

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_velocity_correction.h"
#include "fix_property_atom.h"

#include "atom.h"
#include "compute_pair_gran_local.h"
#include "fix_property_particle.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "mech_param_gran.h"
#include "neigh_list.h"
#include "pair_gran.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingVelocityCorrection::FixCfdCouplingVelocityCorrection(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg)
{       
    if(comm->me==0)printf("Velocity correction for periodic flow is activated\n");
    
    if (strcmp(this->style,"couple/cfd/velocityCorrection") != 0) 
    {
        error->fix_error(FLERR,this,"unknown keyword");
    }

    fix_coupling = NULL;
    fix_velocity_correction= NULL;

    velCorr_total[0] = 0.0;
    velCorr_total[1] = 0.0;
    velCorr_total[2] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingVelocityCorrection::~FixCfdCouplingVelocityCorrection()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingVelocityCorrection::post_create()
{
          
     fix_velocity_correction = static_cast<FixPropertyParticle*>(modify->find_fix_property("velocity_correction","property/particle","vector",3,0,this->style,false));
     
     if(!fix_velocity_correction)
     {
        const char* fixarg[11];
        fixarg[0]="velocity_correction";
        fixarg[1]="all";
        fixarg[2]="property/particle";
        fixarg[3]="velocity_correction";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_velocity_correction = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);	 
     }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingVelocityCorrection::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_velocity_correction) modify->delete_fix("velocity_correction");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingVelocityCorrection::setmask()
{
    int mask = 0;
    //mask |= POST_FORCE | PRE_FORCE;
    mask |= PRE_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingVelocityCorrection::init()
{    
     if(modify->n_fixes_style(style) != 1)
         error->fix_error(FLERR,this,"More than one fix of this style is not allowed");
     
      // find coupling fix
     fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
     if(!fix_coupling) error->fix_error(FLERR,this,"Fix couple/cfd/velocityCorrection needs a fix of type couple/cfd");
     
     // Velocity correction from OF     
     fix_coupling->add_pull_property("velocity_correction","vector-atom");     
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingVelocityCorrection::pre_force(int)
{    
     //Empty
}

/* ---------------------------------------------------------------------- */


void FixCfdCouplingVelocityCorrection::post_force(int)
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **velocity_correction = fix_velocity_correction->array_atom;

  // add velocity correction to velocities
  
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        vectorAdd3D(v[i],velocity_correction[i],v[i]);        
    }
  }
}

/* ---------------------------------------------------------------------- */
double FixCfdCouplingVelocityCorrection::compute_vector(int n)
{
  MPI_Sum_Vector(velCorr_total,3,world);
  return velCorr_total[n];
}

