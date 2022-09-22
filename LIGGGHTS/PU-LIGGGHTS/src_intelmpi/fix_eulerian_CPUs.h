/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

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

FixStyle(cpus,FixEulerianCPUs)

#else

#ifndef LMP_EULERIAN_CPUS_H
#define LMP_EULERIAN_CPUS_H

#include "stdio.h"
#include "fix.h"
#include "vector_liggghts.h"
#include <vector>
//#include "polyhedron_search.h" /* deprecated algorithm (JK) */
#include "ray_search.h"


#define PART_ASSIGNED_SIZE 120

namespace LAMMPS_NS {

class FixEulerianCPUs : public Fix{

 public:

  FixEulerianCPUs(class LAMMPS *, int, char **);
  ~FixEulerianCPUs();
 
  virtual void post_create();
  virtual void init();
 
  int setmask();
  virtual void end_of_step();
  
  void updateCpus();
  
  // inline access functions 
  inline int* localtagCFD()
  { return localtagCFD_; } 
  
   int* libraryTagBuffer_;
   
  inline int* libraryTagBuffer()
  { return libraryTagBuffer_; }
  
 
  
  bool debug_flag;
  //list particle global ids (needs to be structured in the same way as the from pointer
  std::vector<int> cpu_indexes;	// DEM to CFD cpu information	
   
  int n_cfd; 
  int n_dem;
  
  inline RaySearch* getRaySearch() const
  {
     return poly_search;
  }
  
  int coord2bin( const double *x, int ic) const;
  
  inline int prevCPU( int i ) const;
  
 protected :
 
  class FixPropertyParticle* fix_prevCPU;
   
 private:
  void reallocPartAssign(int);
  void bin_atoms(int ic);
   
  void bin_atoms_init(); 
  
  // Number of CPUs
  int ncpus_;
    
  // Coordinates for CPUs
  double * coords_;  
  double * bounds__;
  
  void init_search(int);
  
  int* localtagCFD_;
  bool* part_assigned;
  int spart_assigned;
  
  bool is_digit( char c );
  double parse_double( std::string line, int* i, bool& success );
  int parse_int( std::string line );
  vector<double> parse_vector( std::string line, int i );
  void get_face_vector(  vector<double*>& face, vector<double> parsed );
  
  //PolyhedronSearch* poly_search; /* deprecated algorithm (JK) */
  RaySearch* poly_search;
  
  double tol;
  
  bool compress_flag;
  bool polyhedron_flag;
  
};

}

#endif
#endif
