/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Class to provide scaling due to stiffness effects
   Contributing author:
   Ali Ozel (Princeton University)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "modify.h"
#include "domain.h"
#include "force.h"
#include "neighbor.h"
#include "error.h"
#include "stiffnessScaling.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

StiffnessScaling::StiffnessScaling(LAMMPS *lmp) : Pointers(lmp)
{
}

/* ---------------------------------------------------------------------- */

StiffnessScaling::~StiffnessScaling()
{
}

/* ----------------------------------------------------------------------
   called as stiffnessScaling command in input script
------------------------------------------------------------------------- */

void StiffnessScaling::command(int narg, char **arg)
{
  // error checks
  // ensure is first in input script, so all commands using length scales
  // can use it

  if(neighbor->skin != 0.3 || domain->box_exist)
    error->all(FLERR,"Illegal stiffnessScaling command, must use this command "
                     "first in input script");

  if(modify->nfix > 0)
    error->all(FLERR,"Illegal stiffnessScaling command, must use this command "
                     "before any 'fix' command");

  // parse arguments

  if (narg < 1) error->all(FLERR,"Illegal stiffnessScaling command");

  int iarg = 0;
  double chiStiff = force->numeric(FLERR,arg[iarg++]);
  if(chiStiff < 1.)
    error->all(FLERR,"Illegal stiffnessScaling command, chiStiff > 1 expected");

  //printf(" stiffnessScaling chiStiff %f \n", chiStiff);

  // set stiffnessScaling in force class
  force->chiStiff = chiStiff;

  while(iarg < narg)
  {
      if(strcmp(arg[iarg],"model_check") == 0) {
          if (narg < iarg+2)
            error->all(FLERR,"Illegal stiffnessScaling command, not enough arguments");
          if(strcmp(arg[iarg+1],"error") == 0)
            force->error_stiffScaling = true;
          else if(strcmp(arg[iarg+1],"warn") == 0)
            force->error_stiffScaling = false;
          else
            error->all(FLERR,"Illegal stiffnessScaling command, expecting 'error' or 'warn' after 'model_check'");
          iarg += 2;
      }
      else
        error->all(FLERR,"Illegal stiffnessScaling command, unknown keyword");
  }
}
