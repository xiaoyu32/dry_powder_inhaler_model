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
#include <stdlib.h>
#include "math_extra_liggghts.h"
#include "math_extra_dist_lineSegment.h"

namespace MathExtraDist {

// **********************************************************************
lineSegment::lineSegment(double *line, double *cLine, double *segment, double *cSegment, double segmentExtent)
:
 line_(line),
 cLine_(cLine),
 segment_(segment),
 cSegment_(cSegment),
 segmentExtent_(segmentExtent)
{
  mLineParameter_       = -1e10;    
  mSegmentParameter_    = -1e10;   
  mClosestPoint0_[0] = mClosestPoint0_[1] = mClosestPoint0_[2] = -1e10;
  mClosestPoint1_[0] = mClosestPoint1_[1] = mClosestPoint1_[2] = -1e10;

}

// **********************************************************************
lineSegment::~lineSegment(){}

// **********************************************************************
//    MEMBER FUNCTIONS 
// **********************************************************************
double lineSegment::GetSquared()
{

    double diff[3];
    for (int iDir=0; iDir<3; iDir++)
        diff[iDir] = cLine_[iDir] - cSegment_[iDir];

    double a01 = -1.0 * LAMMPS_NS::vectorDot3D( line_, segment_);
    double b0 = LAMMPS_NS::vectorDot3D(diff, line_);
    double c = diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];

    double det = MathExtraLiggghts::abs( 1.0 - a01*a01 );
    double b1, s0, s1, sqrDist, extDet;

    if (det >= SMALL_TRIMESH)
    {
        // The line and segment are not parallel.
        b1 = -1.0 * LAMMPS_NS::vectorDot3D(diff, segment_);
        s1 = a01*b0 - b1;
        extDet = segmentExtent_*det;

        if (s1 >= -extDet)
        {
            if (s1 <= extDet)
            {
                // Two interior points are closest, one on the line and one
                // on the segment.
                double invDet = 1.0/det;
                s0  = (a01*b1 - b0)*invDet;
                s1 *= invDet;
                sqrDist = s0*(s0     + a01*s1 + 2.0*b0)
                        + s1*(a01*s0 + s1     + 2.0*b1)
                        + c;
            }
            else
            {
                // The endpoint e1 of the segment and an interior point of
                // the line are closest.
                s1 = segmentExtent_;
                s0 = -(a01*s1 + b0);
                sqrDist = -s0*s0 
                        + s1*(s1 + 2.0*b1) 
                        + c;
            }
        }
        else
        {
            // The end point e0 of the segment and an interior point of the
            // line are closest.
            s1 = -segmentExtent_;
            s0 = -(a01*s1 + b0);
            sqrDist = -s0*s0
                    + s1*(s1 + 2.0*b1) 
                    + c;
        }
    }
    else
    {
        // The line and segment are parallel.  Choose the closest pair so that
        // one point is at segment center.
        s1 = 0.0;
        s0 = -b0;
        sqrDist = b0*s0 
                + c;
    }

    for (int iDir=0; iDir<3; iDir++)
    {
        mClosestPoint0_[iDir]    = cLine_[iDir]    + s0 * line_[iDir];
        mClosestPoint1_[iDir]    = cSegment_[iDir] + s1 * segment_[iDir];
    }
    mLineParameter_    = s0;
    mSegmentParameter_ = s1;

    // Account for numerical round-off errors.
    if (sqrDist < 0.0)
        sqrDist = 0.0;

/*
    printf("lineSegment::GetSquared mLineParameter_: %g, mSegmentParameter_: %g, mClosestPoint0: %g %g %g, mClosestPoint1: %g %g %g. \n",
                  mLineParameter_, mSegmentParameter_,
                  mClosestPoint0_[0],mClosestPoint0_[1],mClosestPoint0_[2],
                  mClosestPoint1_[0],mClosestPoint1_[1],mClosestPoint1_[2]
              );
*/

    return sqrDist;

}


}
