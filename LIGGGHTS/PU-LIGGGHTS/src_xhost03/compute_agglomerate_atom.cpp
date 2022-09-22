/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include "compute_agglomerate_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "group.h"

using namespace LAMMPS_NS;

enum{AGGLOMERATE,MASK,COORDS};

/* ---------------------------------------------------------------------- */

ComputeAgglomerateAtom::ComputeAgglomerateAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  agglomerateID(NULL)
{
  if (narg != 4) error->all(FLERR,"Illegal compute agglomerate/atom command");

  double cutoff = force->numeric(FLERR,arg[3]);
  cutsq = cutoff*cutoff;

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_forward = 3;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeAgglomerateAtom::~ComputeAgglomerateAtom()
{
  memory->destroy(agglomerateID);
}

/* ---------------------------------------------------------------------- */

void ComputeAgglomerateAtom::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute agglomerate/atom unless atoms have IDs");
  if (force->pair == NULL)
    error->all(FLERR,"Compute agglomerate/atom requires a pair style to be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,
               "Compute agglomerate/atom cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list
  // full required so that pair of atoms on 2 procs both set their agglomerateID

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"agglomerate/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute agglomerate/atom");

}

/* ---------------------------------------------------------------------- */

void ComputeAgglomerateAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeAgglomerateAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,rtmp;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow agglomerateID array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(agglomerateID);
    nmax = atom->nmax;
    memory->create(agglomerateID,nmax,"agglomerate/atom:agglomerateID");
    vector_atom = agglomerateID;
  }

  // invoke full neighbor list (will copy or build if necessary)
  // on the first step of a run, set preflag to one in neighbor->build_one(...)

  if (update->firststep == update->ntimestep) neighbor->build_one(list,1);
  else neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // if update->post_integrate set:
  // a dynamic group in FixGroup is invoking a variable with this compute
  // thus ghost atom coords need to be up-to-date after initial_integrate()

  if (update->post_integrate) {
    commflag = COORDS;
    comm->forward_comm_compute(this);
  }

  // if group is dynamic, insure ghost atom masks are current

  if (group->dynamic[igroup]) {
    commflag = MASK;
    comm->forward_comm_compute(this);
  }

  // every atom starts in its own agglomerate, with agglomerateID = atomID

  tagint *tag = atom->tag;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) agglomerateID[i] = tag[i];
    else agglomerateID[i] = 0;
  }

  // loop until no more changes on any proc:
  // acquire agglomerateIDs of ghost atoms
  // loop over my atoms, checking distance to neighbors
  // if both atoms are in agglomerate, assign lowest agglomerateID to both
  // iterate until no changes in my atoms
  // then check if any proc made changes

  commflag = AGGLOMERATE;
  double **x = atom->x;
  double *r = atom->radius;

  int change,done,anychange;

  while (1) {
    comm->forward_comm_compute(this);

    change = 0;
    while (1) {
      done = 1;
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;

        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
	rtmp = r[i];
	
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          if (!(mask[j] & groupbit)) continue;
          if (agglomerateID[i] == agglomerateID[j]) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
	  

          rsq = delx*delx + dely*dely + delz*delz - (r[j] + rtmp)*(r[j] + rtmp);
	  
	  //printf("rsq is %f \n",rsq);
	  //printf("cutsq is %f \n",cutsq);

          if (rsq < cutsq) {
            agglomerateID[i] = agglomerateID[j] = MIN(agglomerateID[i],agglomerateID[j]);
            done = 0;
          }
        }
      }
      if (!done) change = 1;
      if (done) break;
    }

    // stop if all procs are done

    MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
    if (!anychange) break;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeAgglomerateAtom::pack_forward_comm(int n, int *list, double *buf,
                                          int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (commflag == AGGLOMERATE) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = agglomerateID[j];
    }
  } else if (commflag == MASK) {
    int *mask = atom->mask;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(mask[j]).d;
    }
  } else if (commflag == COORDS) {
    double **x = atom->x;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeAgglomerateAtom::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == AGGLOMERATE) {
    for (i = first; i < last; i++) agglomerateID[i] = buf[m++];
  } else if (commflag == MASK) {
    int *mask = atom->mask;
    for (i = first; i < last; i++) mask[i] = (int) ubuf(buf[m++]).i;
  } else if (commflag == COORDS) {
    double **x = atom->x;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeAgglomerateAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
