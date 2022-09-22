/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2015-     DCS Computing GmbH, Linz
   Copyright 2015-     TU Graz

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
#include "group.h"
#include "comm.h"
#include "math.h"
#include "vector_liggghts.h"
#include "fix_cfd_coupling_convection_impl.h"
#include "fix_property_particle.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvectiveImpl::FixCfdCouplingConvectiveImpl(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
    fix_coupling = NULL;
    fix_heatFluid  = fix_heatTransCoeff = NULL;

    int iarg = 3;

}

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvectiveImpl::~FixCfdCouplingConvectiveImpl()
{
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectiveImpl::pre_delete(bool unfixflag)
{
    if(fix_heatFluid)       modify->delete_fix("fix_heatFluid");
    if(fix_heatTransCoeff)  modify->delete_fix("fix_heatTransCoeff");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingConvectiveImpl::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectiveImpl::post_create()
{
  //  register fluid temperature and transfer coefficient
  if(!fix_heatFluid)
  {
        const char* fixarg[11];
        fixarg[0]="heatFluid";
        fixarg[1]="all";
        fixarg[2]="property/particle";
        fixarg[3]="heatFluid";
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
        fix_heatFluid = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  if(!fix_heatTransCoeff)
  {
        const char* fixarg[11];
        fixarg[0]="heatTransCoeff";
        fixarg[1]="all";
        fixarg[2]="property/particle";
        fixarg[3]="heatTransCoeff";
        fixarg[4]="scalar";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
        fix_heatTransCoeff = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvectiveImpl::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->all(FLERR,"More than one fix of style couple/cfd/convectiveImpl is not allowed");

    // find coupling fix
    fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling)
      error->all(FLERR,"Fix couple/cfd/convectiveImpl needs a fix of type couple/cfd");

    //values to send to OF, this is the SURFACE TEMPERATURE!
    fix_coupling->add_push_property("Temp","scalar-atom");

    //values to come from OF
    fix_coupling->add_pull_property("heatFluid","scalar-atom");
    fix_coupling->add_pull_property("heatTransCoeff","scalar-atom");

    fix_heatFluid = static_cast<FixPropertyParticle*>(modify->find_fix_property("heatFluid","property/particle","scalar",0,0,style));

    fix_heatTransCoeff = static_cast<FixPropertyParticle*>(modify->find_fix_property("heatTransCoeff","property/particle","scalar",0,0,style));
}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingConvectiveImpl::post_force(int)
{

    //This function may only be used for debug!
/*
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *fluidTemp      = fix_heatFluid->vector_atom;
  double *fluidHeatCoeff = fix_heatTransCoeff->vector_atom;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      printf("fluidTemp/HeatCoeff[%d]: %g %g \n",
             i, fluidTemp[i], fluidHeatCoeff[i]
            );
*/
}
