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
   Richard Berger (JKU Linz)
   Philippe Seil (JKU Linz)
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
------------------------------------------------------------------------- */

#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface.h"
#include "fix_property_particle.h"
#include "modify.h"
#include "container.h"
#include "bounding_box.h"
#include "nbin.h"
#include "atom.h"
#include "domain.h"
#include "vector_liggghts.h"
#include "update.h"
#include <stdio.h>
#include <algorithm>
#include "pair_gran.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "force.h"
#include "npair.h"


#include "bin_hash_table.h"


using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL_DELTA skin/(70.*M_PI)

FixNeighlistMesh::FixNeighlistMesh(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp,narg,arg),
  fix_nneighs_(0),
  fix_nneighs_name_(0),
  buildNeighList(false),
  numAllContacts_(0),
  globalNumAllContacts_(false),
  mbinx(0),
  mbiny(0),
  mbinz(0),
  maxhead(0),
  bins(NULL),
  binhead( NULL ), //JK, changed to reference
  skin(0.0),
  distmax(0.0),
  x(NULL),
  r(NULL),
  changingMesh(false),
  changingDomain(false),
  last_bin_update(-1)
{

    pair_gran_ = NULL;

    if(!modify->find_fix_id(arg[3]) || !dynamic_cast<FixMeshSurface*>(modify->find_fix_id(arg[3])))
        error->fix_error(FLERR,this,"illegal caller");

    caller_ = static_cast<FixMeshSurface*>(modify->find_fix_id(arg[3]));
    mesh_ = caller_->triMesh();

    groupbit_wall_mesh = groupbit;
}

/* ---------------------------------------------------------------------- */

FixNeighlistMesh::~FixNeighlistMesh()
{
    delete [] fix_nneighs_name_;
    last_bin_update = -1;
}


void FixNeighlistMesh::init()
{
    if(!force->pair_match("gran", 0))
	error->fix_error(FLERR,this,"Please use a granular pair style for fix contacthistory");
    pair_gran_ = static_cast<PairGran*>(force->pair_match("gran", 0));
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::post_create()
{
    // register
    if(!fix_nneighs_)
    {
        const char* fixarg[9];
        delete [] fix_nneighs_name_;
        fix_nneighs_name_ = new char[strlen(mesh_->mesh_id())+1+14];
        sprintf(fix_nneighs_name_,"n_neighs_mesh_%s",mesh_->mesh_id());

        fixarg[0]=fix_nneighs_name_;
        fixarg[1]="all";
        fixarg[2]="property/particle";
        fixarg[3]=fix_nneighs_name_;
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart - REQUIRED!
        fixarg[6]="yes";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_nneighs_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);

        fix_nneighs_->just_created = false;
    }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::initializeNeighlist()
{
    changingMesh = mesh_->isMoving() || mesh_->isDeforming();
    changingDomain = (domain->nonperiodic == 2) || domain->box_change;

    // remove old lists, init new ones
    
    const size_t nall = mesh_->sizeLocal()+mesh_->sizeGhost();

    while(triangles.size() > nall) {
        triangles.pop_back();
    }

    while(triangles.size() < nall) {
        triangles.push_back(TriangleNeighlist());
    }

    for(size_t iTri = 0; iTri < nall; iTri++) {
        TriangleNeighlist & triangle = triangles[iTri];
        triangle.contacts.reserve(std::max(triangle.contacts.capacity(), static_cast<size_t>(128)));
    }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::setup_pre_force(int foo)
{
    
    pre_neighbor();
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::min_setup_pre_force(int foo)
{
    
    pre_neighbor();
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_delete(bool unfixflag)
{
    if(unfixflag)
    {
        modify->delete_fix(fix_nneighs_->id);
    }
}

/* ---------------------------------------------------------------------- */

int FixNeighlistMesh::setmask()
{
    int mask = 0;
    mask |= MIN_PRE_NEIGHBOR;
    mask |= PRE_NEIGHBOR;
    mask |= MIN_PRE_FORCE;
    mask |= PRE_FORCE;
    mask |= POST_RUN;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_neighbor()
{
    buildNeighList = true;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::min_pre_force(int vflag)
{
    
    pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_force(int)
{
    
    if(!buildNeighList) return;

    changingMesh = mesh_->isMoving() || mesh_->isDeforming();
    changingDomain = (domain->nonperiodic == 2) || domain->box_change;

    buildNeighList = false;
    numAllContacts_ = 0;

    // set num_neigh = 0
    memset(fix_nneighs_->vector_atom, 0, atom->nlocal*sizeof(double));

    x = atom->x;
    r = atom->radius;

    if(neighbor->style != 1)
        error->all(FLERR,"Please use style 'bin' in the 'neighbor' command together with triangular walls");

    double rmax = 0.5*(neighbor->cutneighmax - neighbor->skin);
    double prev_skin = skin;
    double prev_distmax = distmax;

    if(changingMesh)
    {
      skin = neighbor->skin; // (XL)    
      distmax = neighbor->cutneighmax + SMALL_DELTA;
    }
    else
    {
      skin = 0.5*neighbor->skin; // (XL)
      distmax = neighbor->cutneighmax - rmax + SMALL_DELTA;
    }

   /* mbinx = neighbor->mbinx;
    mbiny = neighbor->mbiny;
    mbinz = neighbor->mbinz;
    bins = neighbor->bins;
    binhead = neighbor->binhead;
    maxhead = neighbor->maxhead;*/
    
    //TODO: understand listIndex, bins 
    const int listIndex = pair_gran_->list->index;
    
    mbinx = neighbor->neigh_bin[listIndex]->mbinx;
    mbiny = neighbor->neigh_bin[listIndex]->mbiny;
    mbinz = neighbor->neigh_bin[listIndex]->mbinz;
    
    bins = neighbor->neigh_bin[listIndex]->bins;
    maxhead = neighbor->neigh_bin[listIndex]->maxbin; // TODO: check maxbin vs maxhead
        
    //binhead = neighbor->neigh_bin[listIndex]->binhead;
    binhead = &(neighbor->neigh_bin[listIndex]->binhead);
    
    
    const size_t nall = mesh_->sizeLocal() + mesh_->sizeGhost();
    
    // update cache if necessary
    if (triangles.size() != nall) 
    {
      initializeNeighlist();
    }

    // update precomputed bins if necessary
    if((skin != prev_skin) || (distmax != prev_distmax) || (neighbor->last_setup_bins > last_bin_update)) 
    {
       generate_bin_list(nall);
    }
        
    // manually trigger binning if no pairwise neigh lists exist
    /*if(0 == neighbor->nlist && bins)
        neighbor->neigh_bin[listIndex]->bin_atoms();
    else if(!bins)
        error->one(FLERR,"wrong neighbor setting for fix neighlist/mesh");
    }*/
    
    
    
    if(0 == neighbor->nlist && binhead){
        neighbor->neigh_bin[listIndex]->bin_atoms();
    }else if(!binhead)
        error->one(FLERR,"wrong neighbor setting for fix neighlist/mesh");
   

    int proc_id;
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );     
    //if( proc_id == 0 ) printf( "Mesh neighbor lists created. (%d)\n", nall );
    
    for(size_t iTri = 0; iTri < nall; iTri++) 
    {
      TriangleNeighlist & triangle = triangles[iTri];
      handleTriangle(iTri);
      numAllContacts_ += triangle.contacts.size();
    }   

    if(globalNumAllContacts_) 
    {
      MPI_Sum_Scalar(numAllContacts_,world);
    }
    
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::handleTriangle(int iTri)
{
    TriangleNeighlist & triangle = triangles[iTri];
    std::vector<int> & neighbors = triangle.contacts;
    int & nchecked = triangle.nchecked;
    int *mask = atom->mask;
    int ixMin(0),ixMax(0),iyMin(0),iyMax(0),izMin(0),izMax(0);
    int nlocal = atom->nlocal;
    double contactDistanceFactor = neighbor->contactDistanceFactor; //(XL)

    neighbors.clear();

    nchecked = 0;

    // only do this if I own particles
    if(nlocal)
    {
      if( changingMesh || changingDomain )
      {
        getBinBoundariesForTriangle(iTri,ixMin,ixMax,iyMin,iyMax,izMin,izMax);
    
        for(int ix=ixMin;ix<=ixMax;ix++) {
          for(int iy=iyMin;iy<=iyMax;iy++) {
            for(int iz=izMin;iz<=izMax;iz++) {
	    
              long iBin = long(iz)*long(mbiny*mbinx) + long(iy*mbinx + ix);
	      
              if(iBin < 0 || iBin >= maxhead) continue;

              int iAtom = (*binhead)[iBin];
              
              while(iAtom != -1 && iAtom < nlocal)
              {
                if(! (mask[iAtom] & groupbit_wall_mesh ))
                {
                    if(bins) iAtom = bins[iAtom];
                    else iAtom = -1;
                    continue;
                }
                nchecked++;

                if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom]*contactDistanceFactor : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
                {
                  neighbors.push_back(iAtom);
                  fix_nneighs_->set_vector_atom_int(iAtom, fix_nneighs_->get_vector_atom_int(iAtom)+1); // num_neigh++
                  
                }
		
                if(bins) iAtom = bins[iAtom];
                else iAtom = -1;
              }
            }
          }
        }
      } else {
      
        const std::vector<long> & triangleBins = triangle.bins;
        const int bincount = triangleBins.size();
	
        for(int i = 0; i < bincount; i++) {
          const long iBin = triangleBins[i];
          
          int iAtom = (*binhead)[iBin];
     
          // only handle local atoms and periodic ghosts -- LIGGGHTS	  
          while(iAtom != -1) // (XL)
          {
             // (XL)
	     if((iAtom > nlocal) && (!domain->is_periodic_ghost(iAtom)))
             {
          	if(bins) iAtom = bins[iAtom];
          	else iAtom = -1;

                continue;
     	     }	              

	    if(! (mask[iAtom] & groupbit_wall_mesh))
            {
                if(bins) iAtom = bins[iAtom];
                else iAtom = -1;
                continue;
            }
	    
            nchecked++;

            if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom]*contactDistanceFactor : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
            {
	      //fprintf( screen, "%d iTri == %d \n", update->ntimestep, iTri );
              neighbors.push_back(iAtom);
              fix_nneighs_->set_vector_atom_int(iAtom, fix_nneighs_->get_vector_atom_int(iAtom)+1); // num_neigh++
            }
	    
            if(bins) iAtom = bins[iAtom];
            else iAtom = -1;
	    
          }
        }
	
      }
    }

}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::getBinBoundariesFromBoundingBox(BoundingBox &b,
      int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax)
{
    double delta = distmax;
    double tri_xmin[3] = {b.xLo-delta,b.yLo-delta,b.zLo-delta};
    double tri_xmax[3] = {b.xHi+delta,b.yHi+delta,b.zHi+delta};
    
    const int listIndex = pair_gran_->list->index;
    
    //int binmin = neighbor->coord2bin(tri_xmin,ixMin,iyMin,izMin);
    //neighbor->neigh_bin[0]->coord2bin(tri_xmin,ixMin,iyMin,izMin);
    neighbor->neigh_pair[listIndex]->coord2bin(tri_xmin,ixMin,iyMin,izMin); 
    
    //int binmax= neighbor->coord2bin(tri_xmax,ixMax,iyMax,izMax);
    //neighbor->neigh_bin[0]->coord2bin(tri_xmax,ixMax,iyMax,izMax);
    neighbor->neigh_pair[listIndex]->coord2bin(tri_xmax,ixMax,iyMax,izMax);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::getBinBoundariesForTriangle(int iTri, int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax)
{
  // disable optimization for movingMesh or shrink-wrapped domain
  if(changingMesh || changingDomain) {
    BoundingBox b = mesh_->getElementBoundingBoxOnSubdomain(iTri);
    // extend bbox by cutneighmax and get bin boundaries
    getBinBoundariesFromBoundingBox(b,ixMin,ixMax,iyMin,iyMax,izMin,izMax);
  } else {
    // use cached boundary information
    const BinBoundary & b = triangles[iTri].boundary;
    ixMin = b.xlo;
    ixMax = b.xhi;
    iyMin = b.ylo;
    iyMax = b.yhi;
    izMin = b.zlo;
    izMax = b.zhi;
  }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::post_run()
{
  last_bin_update = -1; // reset binning for possible next run
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::generate_bin_list(size_t nall)
{
  // precompute triangle bin boundaries
  // disable optimization for changing mesh or domain
  if (!(changingMesh || changingDomain)) {
     
    const int listIndex = pair_gran_->list->index;
        
    double dx = neighbor->neigh_bin[listIndex]->binsizex/2.0; // -- removed division by 2.0 -- (JK)
    double dy = neighbor->neigh_bin[listIndex]->binsizey/2.0; // -- removed division by 2.0 -- (JK)
    double dz = neighbor->neigh_bin[listIndex]->binsizez/2.0; // -- removed division by 2.0 -- (JK)
    double maxdiag = sqrt(dx * dx + dy * dy + dz * dz);
    
    //fprintf( screen, "iAtom = %d (%d %d %d)  nall = %d  distmax = %e \n", (*binhead)[5766], mbinx, mbiny, mbinz, nall, distmax); 
    
    for (size_t iTri = 0; iTri < nall; iTri++) {
      
      TriangleNeighlist & triangle = triangles[iTri];        
      std::vector<long> & binlist = triangle.bins;      
      binlist.clear();
  
      BinBoundary& bb = triangle.boundary;
      BoundingBox b = mesh_->getElementBoundingBoxOnSubdomain(iTri);
            
      // extend bbox by cutneighmax and get bin boundaries
      getBinBoundariesFromBoundingBox(b, bb.xlo, bb.xhi, bb.ylo, bb.yhi, bb.zlo, bb.zhi);

      // look at bins and exclude unnecessary ones
      double center[3] = {0,0,0};
      int total = 0;
      
      /*fprintf( screen, "iTri = %d / %d    (%e) \n", iTri+1, nall, b.xHi );
      
      if( b.xHi >= 0.016 ){
       fprintf( screen, "\n" );
       fprintf( screen, "BB Bound Low = [%d %d %d] \n", bb.xlo, bb.ylo, bb.zlo );
       fprintf( screen, "BB Bound high = [%d %d %d] \n", bb.xhi, bb.yhi, bb.zhi );  
       fprintf( screen, "B Bound Low = [%e %e %e] \n", b.xLo, b.yLo, b.zLo );
       fprintf( screen, "B Bound high = [%e %e %e] \n", b.xHi, b.yHi, b.zHi );  
       
       fprintf( screen, "[%e %e %e] [%e %e %e] [%e %e %e]\n", mesh_->node_(iTri)[0][0], mesh_->node_(iTri)[0][1], mesh_->node_(iTri)[0][2], 
       							      mesh_->node_(iTri)[1][0], mesh_->node_(iTri)[1][1], mesh_->node_(iTri)[1][2],
							      mesh_->node_(iTri)[2][0], mesh_->node_(iTri)[2][1], mesh_->node_(iTri)[2][2]
       );  
       
       fprintf( screen, "\n" );
      }*/
      
      for (int ix = bb.xlo; ix <= bb.xhi; ix++) {
        for (int iy = bb.ylo; iy <= bb.yhi; iy++) {
          for (int iz = bb.zlo; iz <= bb.zhi; iz++) {
            
	    long iBin = long(iz) * long(mbiny * mbinx) + long(iy * mbinx) + long(ix);
	    	    	    
            if (iBin < 0 || iBin >= maxhead)
              continue;

            // determine center of bin (ix, iy, iz)
            neighbor->neigh_bin[listIndex]->bin_center(ix, iy, iz, center);
	    	    	    
            if ( mesh_->resolveTriSphereNeighbuild(iTri, maxdiag, center, distmax + skin) )
            {
               binlist.push_back(iBin);
            }
	    
            total++;
          }
        }
      }
      
      //fprintf( screen, "iTri=%d/%d \n", iTri+1, nall );
      
    }

  }

  last_bin_update = update->ntimestep;
  
  
}

int FixNeighlistMesh::getSizeNumContacts()
{
  return mesh_->sizeLocal() + mesh_->sizeGhost();
}
