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

#include "fix_contact_history.h"
#include "npair_half_bin_newtoff_ghost.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfBinNewtoffGhost::NPairHalfBinNewtoffGhost(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction with partial Newton's 3rd law
   include neighbors of ghost atoms, but no "special neighbors" for ghosts
   owned and ghost atoms check own bin and other bins in stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if i owned and j ghost (also stored by proc owning j)
   pair stored once if i,j are both ghost and i < j
------------------------------------------------------------------------- */

void NPairHalfBinNewtoffGhost::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,which,imol,iatom,moltemplate;
  long ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int xbin,ybin,zbin,xbin2,ybin2,zbin2;
  int *neighptr;

  double **x = atom->x;
  double *radius = atom->radius; //(JK)
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  // -- LIGGGHTS : contact history -- (JK)
  int m,d;
  NeighList *listgranhistory;
  int* touchptr = NULL;
  double *shearptr = NULL;
  int *npartner = NULL,**partner = NULL;
  double **contacthistory = NULL;
  int **firsttouch = NULL;
  double **firstshear = NULL;
  MyPage<int> *ipage_touch = NULL;
  MyPage<double> *dpage_shear = NULL;
  int dnum = list->dnum;
  int nn = 0;
  
  // (JK)
  if( dnum > 0 )
  {
      firstshear = list->firstdouble;
      dpage_shear = list->dpage;
      dpage_shear->reset();
  }

  FixContactHistory *fix_history = list->fix_history; // (JK)
  
  // (JK)
  if (fix_history)
  {
    npartner = fix_history->npartner_;
    partner = fix_history->partner_;
    contacthistory = fix_history->contacthistory_;
    listgranhistory = list->listgranhistory;
  }

  for (i = 0; i < nall; i++) {
    n = 0;
    neighptr = ipage->vget();
    
    // (JK)
    if (fix_history)
    {
       nn = 0;
       shearptr = dpage_shear->vget();
    }
     
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
  
    // when i is a ghost atom, must check if stencil bin is out of bounds
    // only store pair if i < j
    // stores own/own pairs only once
    // stores own/ghost pairs with owned atom only, on both procs
    // stores ghost/ghost pairs only once
    // no molecular test when i = ghost atom

    if (i < nlocal) {
      ibin = atom2bin[i];

      for (k = 0; k < nstencil; k++) {
        for (j = (*binhead)[ibin+stencil[k]]; j >= 0; j = bins[j]) {
          if (j <= i) continue;

          jtype = type[j];
          if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;

	  // --liggghts-- (JK)
	  double radsum = radius[i] + radius[j];
          if (rsq <= cutneighsq[itype][jtype]) {
               if( !fix_history )  neighptr[n] = j;
               else
               {
                  if (rsq <= radsum*radsum)
                  {
                      for (m = 0; m < npartner[i]; m++)
                         if (partner[i][m] == tag[j]) break;
                         if (m < npartner[i])
                         {
                            neighptr[n] = 1;
                            for (d = 0; d < dnum; d++) {shearptr[nn++] = contacthistory[i][m*dnum+d];}
                         }
                         else
                         {
                            neighptr[n] = 0;
                            for (d = 0; d < dnum; d++) {shearptr[nn++] = 0.0;}     
                         }
                 }
                 else
                 {
                   neighptr[n] = 0;
                   for (d = 0; d < dnum; d++) {shearptr[nn++] = 0.0;} 
                 }
              }
              n++; // (JK)
	  }
        }
      }

    } else {
      ibin = coord2bin(x[i],xbin,ybin,zbin);
      for (k = 0; k < nstencil; k++) {
        xbin2 = xbin + stencilxyz[k][0];
        ybin2 = ybin + stencilxyz[k][1];
        zbin2 = zbin + stencilxyz[k][2];
        if (xbin2 < 0 || xbin2 >= mbinx ||
            ybin2 < 0 || ybin2 >= mbiny ||
            zbin2 < 0 || zbin2 >= mbinz) continue;
        for (j = (*binhead)[ibin+stencil[k]]; j >= 0; j = bins[j]) {
          if (j <= i) continue;

          jtype = type[j];
          if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;

          // --liggghts-- (JK)
          double radsum = radius[i] + radius[j];
          if (rsq <= cutneighghostsq[itype][jtype])
          {
               if( !fix_history )  {neighptr[n] = j;}
               else
               {
                  if (rsq <= radsum*radsum)
                  {
                      for (m = 0; m < npartner[i]; m++)
                         if (partner[i][m] == tag[j]) break;
                         if (m < npartner[i])
                         {
                            neighptr[n] = 1;
                            for (d = 0; d < dnum; d++) {shearptr[nn++] = contacthistory[i][m*dnum+d];}
                         }
                         else
                         {
                            neighptr[n] = 0;
                            for (d = 0; d < dnum; d++) {shearptr[nn++] = 0.0;}
                         }
                 }
                 else
                 {
                   neighptr[n] = 0;
                   for (d = 0; d < dnum; d++) {shearptr[nn++] = 0.0;}
                 }
              }
              n++; // (JK)
          }

        }
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    if( dnum > 0 )
    {
        firstshear[i] = shearptr;
        dpage_shear->vgot(nn);
    }
  }
 
  list->inum = atom->nlocal;
  list->gnum = inum - atom->nlocal;
}
