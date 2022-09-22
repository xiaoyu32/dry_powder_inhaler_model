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

#ifdef COMPUTE_CLASS

ComputeStyle(agglomerate/atom,ComputeAgglomerateAtom)

#else

#ifndef LMP_COMPUTE_AGGLOMERATE_ATOM_H
#define LMP_COMPUTE_AGGLOMERATE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeAgglomerateAtom : public Compute {
 public:
  ComputeAgglomerateAtom(class LAMMPS *, int, char **);
  ~ComputeAgglomerateAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();

 private:
  int nmax,commflag;
  double cutsq;
  class NeighList *list;
  double *agglomerateID;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use compute agglomerate/atom unless atoms have IDs

Atom IDs are used to identify agglomerates.

E: Compute agglomerate/atom requires a pair style to be defined

This is so that the pair style defines a cutoff distance which
is used to find agglomerates.

E: Compute agglomerate/atom cutoff is longer than pairwise cutoff

Cannot identify agglomerates beyond cutoff.

W: More than one compute agglomerate/atom

It is not efficient to use compute agglomerate/atom  more than once.

*/
