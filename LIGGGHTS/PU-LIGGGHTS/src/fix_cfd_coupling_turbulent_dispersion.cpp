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
//#include <random>
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_turbulent_dispersion.h"
#include "fix_property_particle.h"

#include <cmath>
#include <cfloat>


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingTurbulentDispersion::FixCfdCouplingTurbulentDispersion(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg),
    fix_coupling_(0),
    fix_vgfluc_(0),
    fix_shearRate_(0)
    
{
    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"fluidViscosity") == 0) {
            if(narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'fluidViscosity'");
	
	    this->nu_f = atof(arg[iarg+1]);
	
            iarg++;
            hasargs = true;
        }
	
	if(strcmp(arg[iarg],"delta") == 0) {
            if(narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'delta'");
	
	    this->delta = atof(arg[iarg+1]);
	
            iarg++;
            hasargs = true;
        }
	
	++iarg;

    }

}

/* ---------------------------------------------------------------------- */

FixCfdCouplingTurbulentDispersion::~FixCfdCouplingTurbulentDispersion()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingTurbulentDispersion::post_create()
{
    // register fluctuating gas velocity
    fix_vgfluc_ = static_cast<FixPropertyParticle*>(modify->find_fix_property("vgFluc","property/particle","vector",3,0,this->style,false));
    if (!fix_vgfluc_)
    {
    	const char* fixarg[11];
    	fixarg[0]="vgFluc";
    	fixarg[1]="all";
    	fixarg[2]="property/particle";
    	fixarg[3]="vgFluc";
    	fixarg[4]="vector"; // 1 vector per particle to be registered
    	fixarg[5]="yes";    // restart
    	fixarg[6]="no";     // communicate ghost
    	fixarg[7]="no";     // communicate rev
    	fixarg[8]="0.";
    	fixarg[9]="0.";
    	fixarg[10]="0.";
    	fix_vgfluc_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
    
    // register shear rate
    fix_shearRate_ = static_cast<FixPropertyParticle*>(modify->find_fix_property("shearRate","property/particle","scalar",1,0,this->style,false));
    if (!fix_shearRate_)
    {
    	const char* fixarg[9];
    	fixarg[0]="shearRate";
    	fixarg[1]="all";
    	fixarg[2]="property/particle";
    	fixarg[3]="shearRate";
    	fixarg[4]="scalar"; // 1 vector per particle to be registered
    	fixarg[5]="yes";    // restart
    	fixarg[6]="no";     // communicate ghost
    	fixarg[7]="no";     // communicate rev
    	fixarg[8]="0.";
    	fix_shearRate_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
     }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingTurbulentDispersion::pre_delete(bool unfixflag)
{
    if(unfixflag) modify->delete_fix("vgFluc");
    if(unfixflag) modify->delete_fix("shearRate");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingTurbulentDispersion::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingTurbulentDispersion::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"Fix couple/cfd/force needs a fix of type couple/cfd");
    
    fix_vgfluc_ = static_cast<FixPropertyParticle*>(modify->find_fix_property("vgFluc","property/particle","vector",3,0,this->style,true));
    fix_shearRate_ = static_cast<FixPropertyParticle*>(modify->find_fix_property("shearRate","property/particle","scalar",1,0,this->style,true));

    //  values to be transfered to OF
    fix_coupling_->add_push_property("vgFluc","vector-atom");

    // values to come from OF  
    fix_coupling_->add_pull_property("shearRate","scalar-atom");
  
}

inline double rnorm()
{
	const double U1 = ((double) rand()/(RAND_MAX));
	const double U2 = ((double) rand()/(RAND_MAX));
	
	return std::min(std::max(sqrt(-2.0*log(U1))*cos(2.0*3.1415926*U2),-3.72),3.72);
}

void FixCfdCouplingTurbulentDispersion::post_force(int)
{
  int nlocal = atom->nlocal;
  
  double *S = fix_shearRate_->vector_atom;
  double **vgFluc = fix_vgfluc_->array_atom;
  double nu_e, eps, tau_sgs, k_sgs, sigma, a1, a2;

  const double Cs = 0.08;
  const double CL = 0.4819;
  double dt = update->dt;

  for (int i = 0; i < nlocal; i++)
  {
	nu_e = Cs*delta*Cs*delta*abs(S[i]);
	eps = (nu_f + nu_e)*abs(S[i]);
	k_sgs = 5*Cs*Cs*delta*delta*S[i]*S[i];
	tau_sgs = CL*k_sgs/std::max(eps,1e-10);
	
	a1 = exp(-dt/std::max(tau_sgs,1e-10));
	a2 = sqrt(1 - a1*a1);
	sigma = sqrt(2./3.*k_sgs);

	for (int j = 0; j < 3; j++)
	{
		fix_vgfluc_->array_atom[i][j] = a1*vgFluc[i][j] + a2*sigma*rnorm();
	}
	
  }
  
      
}

