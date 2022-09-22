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
   
   Struct to define a segment (= a spherocylinder, i.e., a cylinder with 
   two hemispherical end caps)
------------------------------------------------------------------------- */
#ifndef SEGMENT_INTERFACE_H_
#define SEGMENT_INTERFACE_H_

#include <string>

namespace LIGGGHTS {
namespace ContactModels {

struct SegmentData {
  double center[3];
  double orientation[3];
  
  double length;    //EXCLUDING the hemispherical end caps
  double radius;    
  
  double segmentParameter;  //indicating the position within the segment of current interest
};

}
}

#endif /* CONTACT_INTERFACE_H_ */
