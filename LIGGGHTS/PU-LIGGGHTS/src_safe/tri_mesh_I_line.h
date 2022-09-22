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
   Contributing authors:
   Stefan Radl, TU Graz
------------------------------------------------------------------------- */

#ifndef LMP_TRI_MESH_I_LINE_H
#define LMP_TRI_MESH_I_LINE_H

#define VERBOSE 0 //developer to activate/de-activate for testing and debugging

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveTriSegmentContact( int iPart,int nTri, 
                                                   double *line, double *cLine, double length, double cylRadius,      //input: line properties
                                                   double *delta, double &segmentParameter    //output
                                                 )
  {
    // this is the overlap algorithm, neighbor list build is separate

    double bary[3];
    return resolveTriSegmentContactBary(iPart,nTri,line,cLine,length,cylRadius,delta,segmentParameter,bary);
  }

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveTriSegmentContactBary( int iPart,     int nTri,                          
                                                       double *line,  double *cLine, double length, double cylRadius,     //input: line properties
                                                       double *delta, double &segmentParameter, double *bary                //output
                                                     )
  {
    double **n = node_(nTri);
    double distMinusRadius;

    //basic temporary variables
    double extent = 0.5 * length;
    double mClosestPoint1[3]; //the coordinates of the closest points
    bool*  edgeActiveList = edgeActive(nTri);

#if VERBOSE
    printf("...resolveTriSegmentContactBary, nTri: %d, nodes: (%g %g %g) (%g %g %g) (%g %g %g) \n", 
              nTri,
              n[0][0],n[0][1],n[0][2],
              n[1][0],n[1][1],n[1][2],
              n[2][0],n[2][1],n[2][2]
          );
#endif

    //generate object for line-triangle calculation
    double *surfNorm = SurfaceMeshBase::surfaceNorm(nTri);
    MathExtraDist::lineTriangle queryLT(line, cLine, extent, n, surfNorm,edgeActiveList);
    double sqrDist            = queryLT.GetSquared();
           segmentParameter   = queryLT.GetLineParameter();
    double pointOnSegment[3];

    if (segmentParameter >= -extent)
    {
        if (segmentParameter <= extent)
        {
            queryLT.GetClosestPoint0(pointOnSegment);  //some point within the line segment
            
            //if line was parallel, handle it as a sphere-wall interaction
            if( sqrDist<(-LARGE_TRIMESH) && MathExtraLiggghts::abs(segmentParameter)<SMALL_TRIMESH )
            {
                distMinusRadius =
                 resolveTriSphereContactBary(iPart,nTri, cylRadius, pointOnSegment,
                                             delta,bary);
            }
            else
            {
                queryLT.GetClosestPoint1(mClosestPoint1);  //point on the triangle
                bary[0] = queryLT.GetTriangleBary(0);
                bary[1] = queryLT.GetTriangleBary(1);
                bary[2] = queryLT.GetTriangleBary(2);
                distMinusRadius = sqrt(sqrDist) - cylRadius;
                calcDist(pointOnSegment,mClosestPoint1,delta);
            }
        }
        else
        {
            segmentParameter = extent;
            for (int iDir=0; iDir<3; iDir++)
                pointOnSegment[iDir] = cLine[iDir] + extent*line[iDir]; //end of line segment point

            //Easy - Use algorithm for sphere contact
            distMinusRadius =
              resolveTriSphereContactBary(iPart,nTri, cylRadius, pointOnSegment,
                                          delta,bary);
        }
    }
    else
    {
            segmentParameter = -extent;
            for (int iDir=0; iDir<3; iDir++)
                pointOnSegment[iDir] = cLine[iDir] - extent*line[iDir]; //start of line segment point 

            //Easy - Use algorithm for sphere contact
            distMinusRadius =
              resolveTriSphereContactBary(iPart,nTri, cylRadius, pointOnSegment,
                                          delta,bary);
    }

#if VERBOSE 
    printf("segmentParameter %g, extent: %g, pointOnSegment: %g %g %g, delta: %g %g %g distMinusRadius: %g. \n", 
            segmentParameter,extent,
            pointOnSegment[0], pointOnSegment[1], pointOnSegment[2],
            delta[0],delta[1],delta[2],
            distMinusRadius
          );
#endif

    return distMinusRadius;
  }

/* ---------------------------------------------------------------------- */

inline bool TriMesh::resolveTriSegmentNeighbuild(int nTri, double *line, double *cLine, double length, double cylRadius,
                                                  double treshold)
{
    //TODO/WARNING: this function is currently NOT called! The user must ensure that 
    //the skin distance is larger enough to build neighbor lists correctly

    double rSphereMax = 0.5 * length + cylRadius; //maximum extension of the spherocylinder from the center
    return resolveTriSphereNeighbuild(nTri, rSphereMax, cLine, treshold);
}


#endif
