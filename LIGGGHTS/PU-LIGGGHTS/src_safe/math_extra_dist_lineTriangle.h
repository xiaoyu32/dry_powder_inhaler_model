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
   Contributing author:
   Stefan Radl, TU Graz

   Parts of this code were derived from the "WildMagic5" package,
   licensed under BOOST:
   http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
------------------------------------------------------------------------- */

#ifndef LMP_MATH_EXTRA_DIST_LINETRIANGLE_H
#define LMP_MATH_EXTRA_DIST_LINETRIANGLE_H

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "error.h"

#define SMALL_TRIMESH 1.e-10
#define LARGE_TRIMESH 1000000

namespace MathExtraDist {

class lineTriangle {

public:
    lineTriangle (double *line, double *cLine, double lineExtend, double **nodes, double *surfNorm, bool *edgeActive);
    ~lineTriangle();

      //Main Function to compute squared distance and closest points
      double GetSquared();

      // Object access
      double GetLineParameter    () const{ return mLineParameter_;    } ;
      double GetSegmentParameter () const{ return mSegmentParameter_; } ;
      int    GetClosestEdgeId    () const{ return closestEdge_;    } ;

      void GetClosestPoint0 (double *output) { output[0] = mClosestPoint0_[0];output[1] = mClosestPoint0_[1];output[2] = mClosestPoint0_[2]; } ;
      void GetClosestPoint1 (double *output) { output[0] = mClosestPoint1_[0];output[1] = mClosestPoint1_[1];output[2] = mClosestPoint1_[2]; } ;

      double GetTriangleBary (int i) const  { return mTriangleBary_[i]; } ;

private:

      double *line_;        //the vector representing the line direction
      double *cLine_;       //the vector representing a point on the line
      double lineExtend_;   //length of the line, only relevant for check with active edges
      double **nodes_;      //the nodes represeting the vertices of the triangle
      double *surfNorm_;    //the surface normal of the triangle
      bool   *edgeActive_;  //list to indicate whether an edge is active or not
      
      int     closestEdge_; //id of the edge which is closest

      // Information about the closest points.
      mutable double mClosestPoint0_[3];
      mutable double mClosestPoint1_[3];

      double mLineParameter_;        
      double mSegmentParameter_;

      double mTriangleBary_[3];     

};


}

#endif
