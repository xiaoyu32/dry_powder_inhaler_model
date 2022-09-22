/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Efield
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

#ifdef FIX_CLASS

FixStyle(print/nlocal,FixNlocal)

#else

#ifndef LMP_FIX_NLOCAL_H
#define LMP_FIX_NLOCAL_H

#include "fix.h"

namespace LAMMPS_NS {
	
   class FixNlocal : public Fix 
   {
       public:	
	   FixNlocal(class LAMMPS*, int, char**);
	   ~FixNlocal();
	   
	   virtual int setmask();
	   virtual void end_of_step();

	   int step;

   };
}
		
#endif
#endif
