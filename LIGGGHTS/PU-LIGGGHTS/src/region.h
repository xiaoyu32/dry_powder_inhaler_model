/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_REGION_H
#define LMP_REGION_H

#include "domain.h"
#include "pointers.h"
#include "math_extra_liggghts.h"

namespace LAMMPS_NS {

class Region : protected Pointers {
 public:
  char *id,*style;
  int interior;                     // 1 for interior, 0 for exterior
  int scaleflag;                    // 1 for lattice, 0 for box
  double xscale,yscale,zscale;      // scale factors for box/lattice units
  double extent_xlo,extent_xhi;     // bounding box on region
  double extent_ylo,extent_yhi;
  double extent_zlo,extent_zhi;
  int bboxflag;                     // 1 if bounding box is computable
  int varshape;                     // 1 if region shape changes over time
  int dynamic;                      // 1 if position/orient changes over time
  int moveflag,rotateflag;          // 1 if position/orientation changes
  int openflag;             // 1 if any face is open
  int open_faces[6];            // flags for which faces are open

  int copymode;                     // 1 if copy of original class

  // contact = particle near region surface (for soft interactions)
  // touch = particle touching region surface (for granular interactions)

  struct Contact {
    double r;                 // distance between particle & surf, r > 0.0
    double delx,dely,delz;    // vector from surface pt to particle
    double radius;            // curvature of region at contact point
    int iwall;            // unique id of wall for storing shear history
    int varflag;              // 1 if wall can be variable-controlled
  };
  Contact *contact;           // list of contacts
  int cmax;                   // max # of contacts possible with region
  int tmax;           // max # of touching contacts possible

  // motion attributes of region
  // public so can be accessed by other classes

  double dx,dy,dz,theta;      // current displacement and orientation
  double v[3];            // translational velocity
  double rpoint[3];       // current origin of rotation axis
  double omega[3];        // angular velocity
  double rprev;               // speed of time-dependent radius, if applicable
  double xcenter[3];          // translated/rotated center of cylinder/sphere (only used if varshape)
  double prev[5];             // stores displacement (X3), angle and if
                              //  necessary, region variable size (e.g. radius)
                              //  at previous time step
  int vel_timestep;           // store timestep at which set_velocity was called
                              //   prevents multiple fix/wall/gran/region calls
  int nregion;                // For union and intersect
  int size_restart;
  int *list;

  Region(class LAMMPS *, int, char **);
  virtual ~Region();
  virtual void init();
  int dynamic_check();

  // called by other classes to check point versus region

  void prematch();
  int match(double, double, double);
  int surface(double, double, double, double);

  virtual void set_velocity();
  void velocity_contact(double *, double *, int);
  virtual void write_restart(FILE *);
  virtual int restart(char *, int&);
  virtual void length_restart_string(int&);
  virtual void reset_vel();

  // implemented by each region, not called by other classes

  virtual int inside(double, double, double) = 0;
  virtual int surface_interior(double *, double) = 0;
  virtual int surface_exterior(double *, double) = 0;
  virtual void shape_update() {}
  virtual void pretransform();
  virtual void set_velocity_shape() {}
  virtual void velocity_contact_shape(double*, double*) {}

  // -- LIGGGHTS METHODS -- (JK)

  // generates a random point within the region
  virtual void generate_random(double *,bool subdomain_flag);

  // reset random gen - is called out of restart by fix that uses region
  void reset_random(int);
  
  // volume calculation based on MC
  virtual void volume_mc(int n_test,bool cutflag,double cut,double &vol_global,double &vol_local);

  // flag if region bbox extends outside simulation domain
  virtual int bbox_extends_outside_box();

  // inside region AND within a minimum distance from surface
  int match_cut(double *,double);

  // inside region OR within a minimum distance from surface
  int match_expandby_cut(double *,double);

  // generate a point inside region AND further away from surface than cut
  virtual void generate_random_shrinkby_cut(double *,double,bool subdomain_flag);

  inline void rand_bounds(bool subdomain_flag, double *lo, double *hi)
  {
      if(!bboxflag) error->one(FLERR,"Impossible to generate random points on region with incomputable bounding box");
      if(subdomain_flag)
      {
          lo[0] = MathExtraLiggghts::max(extent_xlo,domain->sublo[0]);
          lo[1] = MathExtraLiggghts::max(extent_ylo,domain->sublo[1]);
          lo[2] = MathExtraLiggghts::max(extent_zlo,domain->sublo[2]);
          hi[0] = MathExtraLiggghts::min(extent_xhi,domain->subhi[0]);
          hi[1] = MathExtraLiggghts::min(extent_yhi,domain->subhi[1]);
          hi[2] = MathExtraLiggghts::min(extent_zhi,domain->subhi[2]);
          if(lo[0] >= hi[0] || lo[1] >= hi[1] ||lo[2] >= hi[2])
              error->one(FLERR,"Impossible to generate random points on wrong sub-domain");
      }
      else
      {
          vectorConstruct3D(lo,  extent_xlo,extent_ylo,extent_zlo );
          vectorConstruct3D(hi,  extent_xhi,extent_yhi,extent_zhi );
      }

  }


 protected:
  void add_contact(int, double *, double, double, double);
  void options(int, char **);
  void point_on_line_segment(double *, double *, double *, double *);
  void forward_transform(double &, double &, double &);
  double point[3],runit[3];
  
  // -- liggghts --
  int seed;
  class RanPark *random;
  
  
 private:
  char *xstr,*ystr,*zstr,*tstr;
  int xvar,yvar,zvar,tvar;
  double axis[3];

  void inverse_transform(double &, double &, double &);
  void rotate(double &, double &, double &, double);
};

}

#endif

/* ERROR/WARNING messages:

E: Variable name for region does not exist

Self-explanatory.

E: Variable for region is invalid style

Only equal-style variables can be used.

E: Variable for region is not equal style

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region union or intersect cannot be dynamic

The sub-regions can be dynamic, but not the combined region.

E: Region cannot have 0 length rotation vector

Self-explanatory.

*/
