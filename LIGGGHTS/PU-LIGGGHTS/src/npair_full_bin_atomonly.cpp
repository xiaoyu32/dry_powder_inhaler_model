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
#include "npair_full_bin_atomonly.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairFullBinAtomonly::NPairFullBinAtomonly(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void NPairFullBinAtomonly::build(NeighList *list)
{
  int i,j,k,n,itype,jtype;
  long ibin;
  
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag; //(JK)
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
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
  
  
  if( dnum > 0 )
  {
      firstshear = list->firstdouble;
      dpage_shear = list->dpage;
      dpage_shear->reset();
  }
  
  FixContactHistory *fix_history = list->fix_history; 
  
  if (fix_history) 
  {
    npartner = fix_history->npartner_; 
    partner = fix_history->partner_; 
    contacthistory = fix_history->contacthistory_; 
    listgranhistory = list->listgranhistory;
    
    /*firsttouch = listgranhistory->firstneigh;
    firstshear = listgranhistory->firstdouble;
    ipage_touch = listgranhistory->ipage;
    dpage_shear = listgranhistory->dpage;
    dnum = listgranhistory->dnum; */
  }  
  
  //if (fix_history) 
  //{
    //ipage_touch->reset();
    //dpage_shear->reset();
  //}  
  
  //if( dnum > 0 ) printf( "Building granular history list !!!!!!!!!! %d \n", fix_history != NULL  );
  //else printf( "Building list !!!!!!!!!! \n" );
    
  // -- LIGGGHTS : contact history -- (JK)
  //printf( "Here: dnum = %d    list->dnum = %d \n", dnum, list->dnum ); //fixme
  
  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
        
    if (fix_history) 
    {
       nn = 0;
       //touchptr = ipage_touch->vget();
       shearptr = dpage_shear->vget();
    }    
        
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in surrounding bins in stencil including self
    // skip i = j

    ibin = atom2bin[i];

    for (k = 0; k < nstencil; k++) {
       for (j = (*binhead)[ibin+stencil[k]]; j >= 0; j = bins[j]) 
       {
          if (i == j) continue;

          jtype = type[j];
          if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
	  double radsum = radius[i] + radius[j]; // -- liggghts --

          if (rsq <= cutneighsq[itype][jtype]){
	      
	      if( !fix_history )
	      {
	         neighptr[n] = j;
	      }else
	      {
	         //neighptr[n] = j;
		  
        	 if (rsq <= radsum*radsum)
        	 {
		 
        	   for (m = 0; m < npartner[i]; m++)
        	     if (partner[i][m] == tag[j]) break;
        	   if (m < npartner[i]) 
		   {
        	     neighptr[n] = 1;
        	     for (d = 0; d < dnum; d++) 
		     { 
                       shearptr[nn++] = contacthistory[i][m*dnum+d];
        	     }
        	   } else {
        	      neighptr[n] = 0;
        	      for (d = 0; d < dnum; d++) 
		      { 
                	shearptr[nn++] = 0.0;
        	      }
        	   }
		   
        	 } else {
		 
        	   neighptr[n] = 0;
        	   for (d = 0; d < dnum; d++) 
		   { 
                       shearptr[nn++] = 0.0;
        	   }
		   
        	 }
		 
              }

	      // -- LIGGGHTS -- 	    
	      n++;	
	      
	      //printf( "neighptr[n] = %d \n",  neighptr[n] );
	      
	  }
       }
      
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    
    //if( dnum > 0 ) firstdouble[i] = shearptr;
        
    numneigh[i] = n;
    ipage->vgot(n);
    
    if (ipage->status())
    {
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }
    
    if( dnum > 0 )
    {
        firstshear[i] = shearptr;
        dpage_shear->vgot(nn);
    }
    
    /*
    if( fix_history ) 
    {
      //firsttouch[i] = touchptr;
      //ipage_touch->vgot(n);
      //firsttouch[i] = neighptr;
      
      firstshear[i] = shearptr;
      dpage_shear->vgot(nn);
      
      if (ipage_touch->status()) error->one(FLERR,"Neighbor list overflow, boost neigh_modify one (2)");
      if (dpage_shear->status()) error->one(FLERR,"Neighbor list overflow, boost neigh_modify one (3)");
    }    
    */    
  }

  list->inum = inum;
  list->gnum = 0;
  
}
