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
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_TRI_MESH_H
#define LMP_TRI_MESH_H


#define SMALL_TRIMESH 1.e-10
#define LARGE_TRIMESH 1000000

#include "surface_mesh.h"
#include "atom.h"
#include "math_extra_liggghts.h"
#include "math_extra_dist_lineTriangle.h"
#include <fstream>

namespace LAMMPS_NS
{
  
  typedef SurfaceMesh<3,5> SurfaceMeshBase;

  class TriMesh : public SurfaceMeshBase
  {
      public:

        TriMesh(LAMMPS *lmp);
        virtual ~TriMesh();

        double resolveTriSphereContact(int iPart, int nTri, double rSphere, double *cSphere, double *delta);
        double resolveTriSphereContactBary(int iPart, int nTri, double rSphere, double *cSphere,
                                           double *contactPoint,double *bary);
        bool resolveTriSphereNeighbuild(int nTri, double rSphere, double *cSphere, double treshold);

        //Extra for Line Contact Calculation ********
        double resolveTriSegmentContact    (int iPart, int nTri, double *line, double *cLine, double length, double cylRadius,
                                            double *delta, double &segmentParameter);
        double resolveTriSegmentContactBary(int iPart, int nTri, double *line, double *cLine, double length, double cylRadius,
                                            double *delta, double  &segmentParameter, double *bary);
        bool resolveTriSegmentNeighbuild(int nTri, double *line, double *cLine, double length, double cylRadius, double treshold);
        //Extra for Line Contact Calculation ********

        int generateRandomOwnedGhost(double *pos);
        int generateRandomSubbox(double *pos);


        virtual int generateRandomOwnedGhostWithin(double *pos,double delta);
        double calcArea(int n);
	
      protected:

        bool isInElement(double *pos,int i);

      private:

        double calcDist(double *cs, double *closestPoint, double *en0);
        double calcDistToPlane(double *p, double *pPlane, double *nPlane);

        double resolveCornerContactBary(int iTri, int iNode, bool obtuse,
                                    double *p, double *delta, double *bary);
        double resolveEdgeContactBary(int iTri, int iEdge, double *p, double *delta, double *bary);
        double resolveFaceContactBary(int iTri, double *p, double *node0ToSphereCenter, double *delta);


  };

  // *************************************
  #include "tri_mesh_I.h"
  #include "tri_mesh_I_line.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* TRIMESH_H_ */
