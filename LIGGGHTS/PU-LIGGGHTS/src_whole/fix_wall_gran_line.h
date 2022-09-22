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

#ifdef FIX_CLASS
FixStyle(wall/gran/line,FixWallGranLine)
#else

#ifndef LMP_FIX_WALL_GRAN_LINE_H
#define LMP_FIX_WALL_GRAN_LINE_H

#include "fix_wall_gran.h"
#include "fix_mesh_surface.h"
#include "contact_interface.h"
#include "granular_wall.h"
#include <string>
#include <vector>
#include "fix_contact_property_atom_wall.h"
#include "compute_pair_gran_local.h"
#include "segment_interface.h"

namespace LCM = LIGGGHTS::ContactModels;

namespace LAMMPS_NS {

class FixWallGranLine : public FixWallGran {
 public:
  FixWallGranLine(class LAMMPS *, int, char **);
  ~FixWallGranLine();


  void init();

 protected:

   double const* const* kn_;
   double const* const* gamman_;
   double const* const* gammat_;
   double const* const* coefficientFriction_;
   
   class   FixPropertyGlobal* roughnessFix_;
   class   FixPropertyGlobal* stiffnessFactorFix_;
   class   FixPropertyGlobal* fluidViscosityFix_;

   double roughness_;
   double stiffnessFactor_;
   double fluidViscosity_;

   void post_force_mesh(int);

   bool haveLineData_;

   class AtomVecEllipsoid *avec;

   class FixPropertyParticle* fix_orientation_; 

   inline void post_force_eval_contact( LCM::CollisionData & cdata, LCM::SegmentData & mySegment,
                                        double * v_wall, 
                                        int iMesh = -1, FixMeshSurface *fix_mesh = 0, 
                                        TriMesh *mesh = 0, int iTri = 0);
  
   void compute_force(LCM::CollisionData & cdata, LCM::SegmentData & mySegment, double *vwall);

};

}

#endif
#endif
