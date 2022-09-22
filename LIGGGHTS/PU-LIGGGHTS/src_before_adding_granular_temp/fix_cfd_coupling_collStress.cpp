// Collisional stresses fix class 

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
#include "fix_cfd_coupling_collStress.h"
#include "fix_property_particle.h"
#include "compute_pair_gran_local.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "mech_param_gran.h"
#include "neigh_list.h"
#include "pair_gran.h"
#include "compute_stress_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingCollStress::FixCfdCouplingCollStress(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg)
{       
    if(comm->me==0)printf("Collisional stresses will be sent to CFD solver\n");
    
    if (strcmp(this->style,"couple/cfd/collStress") != 0) 
    {
        error->fix_error(FLERR,this,"unknown keyword");
    }

    fix_coupling = NULL;
    fix_collStress= NULL;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingCollStress::~FixCfdCouplingCollStress()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingCollStress::post_create()
{
     
     if(!fix_collStress)
     {
        const char* fixarg[14];
	fixarg[0]="collStress";
        fixarg[1]="all";
        fixarg[2]="property/particle";
        fixarg[3]="collStress";
        fixarg[4]="vector";     //"vector";
        fixarg[5]="no";        // restart
        fixarg[6]="yes";	// communicate ghost
        fixarg[7]="no";         // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fixarg[11]="0.";
        fixarg[12]="0.";
        fixarg[13]="0.";			
	fix_collStress = modify->add_fix_property_atom(14,const_cast<char**>(fixarg),style);		 
     }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingCollStress::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_collStress) modify->delete_fix("collStress");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingCollStress::setmask()
{
    int mask = 0;
    mask = POST_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingCollStress::init()
{    
     if(modify->n_fixes_style(style) != 1)
         error->fix_error(FLERR,this,"More than one fix of this style is not allowed");
     
      // find coupling fix
     fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
     if(!fix_coupling) error->fix_error(FLERR,this,"Fix couple/cfd/collStress needs a fix of type couple/cfd");

      // find collStress
     fix_collStress = static_cast<FixPropertyParticle*>(modify->find_fix_property("collStress","property/particle","vector",0,0,style));
     
     //  values to be transfered to OF
     fix_coupling->add_push_property("collStress","vector-atom");    
     
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingCollStress::pre_force(int)
{    
     //Empty 
}

/* ---------------------------------------------------------------------- */


void FixCfdCouplingCollStress::post_force(int)
{
       
     ComputeStressAtom* compute_stress = static_cast<ComputeStressAtom*> (modify->find_compute_style_strict( "stress/atom", 0 ));
     
     if( !compute_stress ) return;
     
     compute_stress->compute_peratom();
          
     double** stress_ = compute_stress->array_atom;
     double** stress  = fix_collStress->array_atom;
     
     //shallow copy stresses
     
     for( int i = 0; i < atom->nlocal; ++i )
     {
        for( int j = 0; j < 6; ++j )
	{
	    stress[i][j] = stress_[i][j];
	}
     }
     
}


