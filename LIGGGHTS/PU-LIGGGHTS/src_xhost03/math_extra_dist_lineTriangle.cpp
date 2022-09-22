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
#include "math_extra_dist_lineTriangle.h"
#include "math_extra_dist_lineSegment.h"

#define VERBOSE 0 //developer to activate/de-activate for testing and debugging

namespace MathExtraDist {

// **********************************************************************
lineTriangle::lineTriangle(double *line, double *cLine, double lineExtend, double **nodes, double *surfNorm, bool *edgeActive)
:
 line_(line),
 cLine_(cLine),
 lineExtend_(lineExtend),
 nodes_(nodes),
 surfNorm_(surfNorm),
 edgeActive_(edgeActive)
{
  mLineParameter_       = -LARGE_TRIMESH;    
  mSegmentParameter_    = -LARGE_TRIMESH;   
  mClosestPoint0_[0] = mClosestPoint0_[1] = mClosestPoint0_[2] = -LARGE_TRIMESH;
  mClosestPoint1_[0] = mClosestPoint1_[1] = mClosestPoint1_[2] = -LARGE_TRIMESH;

  mTriangleBary_[0] = mTriangleBary_[1] = mTriangleBary_[2] = -LARGE_TRIMESH;

  closestEdge_ = -1;
}

// **********************************************************************
lineTriangle::~lineTriangle()
{

}

// **********************************************************************
//    MEMBER FUNCTIONS 
// **********************************************************************
double lineTriangle::GetSquared()
{ 
    closestEdge_ = -1; //reset
    //Dot the triangle's normal with the line orientation
    double NdotD = LAMMPS_NS::vectorDot3D(surfNorm_, line_); 

    if ( MathExtraLiggghts::abs(NdotD) > SMALL_TRIMESH ) //check if parallel
    {
        // The line and triangle are not parallel, so the line intersects
        // the plane of the triangle.
#if VERBOSE
        printf("lineTriangle: non-parallel line detected, surfNorm: %g %g %g, center: %g %g %g, line: %g %g %g, NdotD: %g \n",
                  surfNorm_[0],surfNorm_[1],surfNorm_[2],
                  cLine_[0],cLine_[1],cLine_[2],
                  line_[0],line_[1],line_[2],
                  NdotD
              );
#endif
        //Compute difference to first node, as well as edges
        double diff[3], edge0[3], edge1[3];
        for (int iDir=0; iDir<3; iDir++)
        {
            diff[iDir]      = cLine_[iDir]    - nodes_[0][iDir];
            edge0[iDir]     = nodes_[1][iDir] - nodes_[0][iDir];
            edge1[iDir]     = nodes_[2][iDir] - nodes_[0][iDir];
        }

        //Generate complement basis
        double u[3],v[3];
        MathExtraLiggghts::generateComplementBasis(u, v, line_);

        //Generate Dot products
        double UdE0   = LAMMPS_NS::vectorDot3D(u,edge0);
        double UdE1   = LAMMPS_NS::vectorDot3D(u,edge1);
        double UdDiff = LAMMPS_NS::vectorDot3D(u,diff);
        double VdE0   = LAMMPS_NS::vectorDot3D(v,edge0);
        double VdE1   = LAMMPS_NS::vectorDot3D(v,edge1);
        double VdDiff = LAMMPS_NS::vectorDot3D(v,diff);
        double invDet = 1.0/(UdE0*VdE1 - UdE1*VdE0);

        // Barycentric coordinates for the point of intersection.
        double b1 = (VdE1*UdDiff - UdE1*VdDiff)*invDet;
        double b2 = (UdE0*VdDiff - VdE0*UdDiff)*invDet;
        double b0 = 1.0 - b1 - b2;

        if ( b0 >= 0.0 && b1 >= 0.0 && b2 >= 0.0 )
        {
            // Line parameter for the point of intersection.
            double DdE0     = LAMMPS_NS::vectorDot3D(line_,edge0);
            double DdE1     = LAMMPS_NS::vectorDot3D(line_,edge1);
            double DdDiff   = LAMMPS_NS::vectorDot3D(line_,diff);
            mLineParameter_ = b1*DdE0 
                            + b2*DdE1 
                            - DdDiff;

            // Barycentric coordinates for the point of intersection.
            mTriangleBary_[0] = b0;
            mTriangleBary_[1] = b1;
            mTriangleBary_[2] = b2;

            // The intersection point is inside or on the triangle.
            for (int iDir=0; iDir<3; iDir++)
            {
                mClosestPoint0_[iDir] = cLine_[iDir] 
                                     + mLineParameter_*line_[iDir];
                mClosestPoint1_[iDir] = nodes_[0][iDir]
                                     + b1*edge0[iDir] 
                                     + b2*edge1[iDir];
            }

            return 0.0;   //the line intersects the triangle, so the distance is zero
        }
    }
    else    //line is exactly parallel. Simply use center position 
            //Thus, return large negative distance to trigger point-wall interaction algorithm
    {
        mLineParameter_ = 0.0;
        for (int iDir=0; iDir<3; iDir++)
            mClosestPoint0_[iDir] = cLine_[iDir];

        return  -2*LARGE_TRIMESH;
        
    }

    // Either (1) the line is not parallel to the triangle and the point of
    // intersection of the line and the plane of the triangle is outside the
    // triangle or (2) the line and triangle are parallel.  Regardless, the
    // closest point on the triangle is on an EDGE of the triangle. Thus, compare
    // the line to all three edges of the triangle.
    double sqrDist = LARGE_TRIMESH;
    int NUM_NODES=3; //Warning: will only work for triangles!
    for(int iEdge=0;iEdge<NUM_NODES;iEdge++)
    {
        int i1 = (iEdge+1)%NUM_NODES; //same numbering as in surface_mesh function 'calcEdgeVecLen'
        int i0 =  iEdge;
       
        //Extract information on the edges of the triangle, and use as "segment"
        double segmentCenter[3];
        double segmentDirection[3];
        double segmentExtent;

        for (int iDir=0; iDir<3; iDir++)
        {
          segmentCenter[iDir]  = 0.5*(  nodes_[i0][iDir]
                                      + nodes_[i1][iDir]
                                     );
          segmentDirection[iDir] =  nodes_[i1][iDir]
                                 -  nodes_[i0][iDir];
        }
        segmentExtent  =    sqrt(   segmentDirection[0]*segmentDirection[0]
                                  + segmentDirection[1]*segmentDirection[1]
                                  + segmentDirection[2]*segmentDirection[2]
                                 ); //this is the full length, divide by 2 later

        //Normalize direction
        for (int iDir=0; iDir<3; iDir++)
            segmentDirection[iDir] /= segmentExtent;

        segmentExtent *= 0.5; //the extend is defined as 1/2 of the size

        //Compute distance of line to edge (i.e., the "segment")
        MathExtraDist::lineSegment queryLS(line_, cLine_, segmentDirection, segmentCenter, segmentExtent);
        double sqrDistTmp       = queryLS.GetSquared();
        double tmpLineParameter = queryLS.GetLineParameter();

        //we must now skip INACTIVE edges in certain situations to avoid double-edge interaction
        //Skip edges that are NOT active and we expect a line-segment collision for this edge
        //this is because we must consider situations in which the
        //endpoint is closest to the edge (i.e., we cannot skip ALL INACTIVE edges)
        if(       !edgeActive_[iEdge] 
             && ( MathExtraLiggghts::abs(tmpLineParameter) <  lineExtend_)
          )
            continue;


        if (sqrDistTmp < sqrDist)
        {
            queryLS.GetClosestPoint0(mClosestPoint0_);  
            queryLS.GetClosestPoint1(mClosestPoint1_);  

            mLineParameter_ = tmpLineParameter;
            sqrDist         = sqrDistTmp;
            
            double    ratio = queryLS.GetSegmentParameter()/segmentExtent;
            mTriangleBary_[i0] = 0.5*(1.0 - ratio);
            mTriangleBary_[i1] = 1.0 - mTriangleBary_[i0];
            mTriangleBary_[3-i0-i1] = 0.0;
            
            closestEdge_ = iEdge;
#if VERBOSE
        printf("lineTriangle: checked closest ACTIVE edge with center %g %g %g and direction: %g %g %g; distTmp: %g, mClosestPoint0: %g %g %g, mLineParameter: %g. \n",
                  segmentCenter[0],segmentCenter[1],segmentCenter[2],
                  segmentDirection[0],segmentDirection[1],segmentDirection[2],
                  sqrt(sqrDistTmp),
                  mClosestPoint0_[0],mClosestPoint0_[1],mClosestPoint0_[2],
                  mLineParameter_
              );
#endif
        }

    }
    return sqrDist;
}



}


