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
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_wall_gran_line.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair_gran.h"
#include "fix_rigid.h"
#include "fix_mesh.h"
#include "fix_contact_history_mesh.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "fix_property_particle.h"
#include "fix_contact_property_atom_wall.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "compute_pair_gran_local.h"
#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface_stress.h"
#include "tri_mesh.h"
#include "primitive_wall.h"
#include "primitive_wall_definitions.h"
#include "mpi_liggghts.h"
#include "neighbor.h"
#include "contact_interface.h"
#include "fix_property_global.h"
#include <vector>
#include "granular_wall.h"
#include "atom_vec_ellipsoid.h"
#include "global_properties.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LAMMPS_NS::PRIMITIVE_WALL_DEFINITIONS;
using namespace LIGGGHTS::Walls;
using namespace LIGGGHTS::ContactModels;

const double SMALL = 1e-12;

#define VERBOSEFixWallGranLine 0 //developer to activate/de-activate for testing and debugging

  // modes for conduction contact area calaculation
  // same as in fix_heat_gran_conduction.cpp
/* ---------------------------------------------------------------------- */

FixWallGranLine::FixWallGranLine(LAMMPS *lmp, int narg, char **arg) :
  FixWallGran(lmp, narg, arg),
  kn_(NULL),
  gamman_(NULL),
  gammat_(NULL),
  coefficientFriction_(NULL),
  roughnessFix_(NULL),
  stiffnessFactorFix_(NULL),
  fluidViscosityFix_(NULL),
  haveLineData_(false),
  fix_orientation_(0)
{

}

/* ---------------------------------------------------------------------- */
FixWallGranLine::~FixWallGranLine()
{
}


/* ---------------------------------------------------------------------- */
void FixWallGranLine::init()
{

     avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
     if (!avec) error->all(FLERR,"FixWallGranLine requires atom style ellipsoid");

     fix_orientation_ = static_cast<FixPropertyParticle*>(modify->find_fix_property("ex","property/particle","vector",0,0,style));
     if(!fix_orientation_)
        error->fix_error(FLERR,this,"Fix FixWallGranLine could NOT find fix with id 'ex'. This is fatal."); 

     //Global Parameters        
     kn_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("kn","property/global","peratomtypepair",atom->ntypes,atom->ntypes,style))->get_array();
     if(kn_==NULL)
        error->fix_error(FLERR,this,"Fix FixWallGranLine could NOT find fix with id 'kn'. This is fatal."); 
        
  gamman_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("gamman","property/global","peratomtypepair",atom->ntypes,atom->ntypes,style))->get_array();
     if(gamman_==NULL)
        error->fix_error(FLERR,this,"Fix FixWallGranLine could NOT find fix with id 'gamman'. This is fatal."); 

  gammat_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("gammat","property/global","peratomtypepair",atom->ntypes,atom->ntypes,style))->get_array();
     if(gammat_==NULL)
        error->fix_error(FLERR,this,"Fix FixWallGranLine could NOT find fix with id 'gammat'. This is fatal."); 
        
  coefficientFriction_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",atom->ntypes,atom->ntypes,style))->get_array();
     if(coefficientFriction_==NULL)
        error->fix_error(FLERR,this,"Fix FixWallGranLine could NOT find fix with id 'coefficientFriction'. This is fatal."); 
  
    roughnessFix_=static_cast<FixPropertyGlobal*>(modify->find_fix_property("roughness","property/global","scalar",0,0,force->pair_style));
    roughness_ = roughnessFix_->compute_scalar(); // const for all types
    if(roughness_<=0.0)
        error->fix_error(FLERR,this,"property/global 'roughness' cannot be smaller or equal 0.0"); 


    stiffnessFactorFix_=static_cast<FixPropertyGlobal*>(modify->find_fix_property("stiffnessFactor","property/global","scalar",0,0,force->pair_style));
    stiffnessFactor_ = stiffnessFactorFix_->compute_scalar(); // const for all types

    if(stiffnessFactor_>1.0)
        error->fix_error(FLERR,this,"property/global 'stiffnessFactor' cannot be larger than 1.0."); 
        

    fluidViscosityFix_=static_cast<FixPropertyGlobal*>(modify->find_fix_property("fluidViscosity","property/global","scalar",0,0,force->pair_style));
    fluidViscosity_ = fluidViscosityFix_->compute_scalar(); // const for all types
  
}


/* ----------------------------------------------------------------------
   post_force for mesh wall
------------------------------------------------------------------------- */

void FixWallGranLine::post_force_mesh(int vflag)
{

    // contact properties
    double v_wall[3],bary[3];
    double delta[3],deltan;
    MultiVectorContainer<double,3,3> *vMeshC;
    double ***vMesh;
    double **orientation = NULL;
    int nlocal = atom->nlocal;
    int nTriAll;

    AtomVecEllipsoid::Bonus *bonus = avec->bonus;
    int     *ellipsoid             = atom->ellipsoid;
    double  *shape;


    orientation = fix_orientation_->array_atom;

    CollisionData cdata;
    cdata.is_wall = true;

    for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
    {

      TriMesh *mesh = FixMesh_list_[iMesh]->triMesh();
      nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
      FixContactHistoryMesh *fix_contact = FixMesh_list_[iMesh]->contactHistory();

      // mark all contacts for delettion at this point
      
      if(fix_contact) fix_contact->markAllContacts();

      if(store_force_contact_)
        fix_wallforce_contact_ = FixMesh_list_[iMesh]->meshforceContact();

      // get neighborList and numNeigh
      FixNeighlistMesh * meshNeighlist = FixMesh_list_[iMesh]->meshNeighlist();

      vectorZeroize3D(v_wall);
      vMeshC = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");

      atom_type_wall_ = FixMesh_list_[iMesh]->atomTypeWall();

      //Main routine for contact and force calculation
      {
        // Compute velocity/rotation tensor of mesh
        if(vMeshC)
           vMesh = vMeshC->begin();


        // loop owned and ghost particles
        for(int iTri = 0; iTri < nTriAll; iTri++)
        {
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();

          for(int iCont = 0; iCont < numneigh; iCont++)
          {
            const int iPart = neighborList[iCont];

            // no need to handle ghost particles
            if(iPart >= nlocal) continue;

            SegmentData mySegment;
            for(int iDir=0;iDir<3;iDir++)
            {
                mySegment.center[iDir]      = x_[iPart][iDir];
                mySegment.orientation[iDir] = orientation[iPart][iDir];
            }
            
            //Extract cylinder parameters
            double major, minor;
            shape   = bonus[ellipsoid[iPart]].shape;
            major   = MathExtraLiggghts::max(shape[0],shape[1],shape[2]);
            minor   = MathExtraLiggghts::min(shape[0],shape[1],shape[2]);
            mySegment.radius = 1.24 * minor / sqrt( log(major/minor) );
            mySegment.length = 2.0*(  major
                                    - mySegment.radius
                                   );  //length EXCLUDING the spherical end regions

            deltan = mesh->resolveTriSegmentContactBary(iPart,iTri, 
                                                    mySegment.orientation, mySegment.center, 
                                                    mySegment.length, mySegment.radius, 
                                                    delta, mySegment.segmentParameter, bary
                                                   );
#if VERBOSEFixWallGranLine
    #include "fix_wall_gran_line_debug.h"
#endif
            //Force calculation
            if(deltan > skinDistance_) //allow force calculation away from the wall
            {
            }
            else
            {
              //WARNING: no contact_history handling for lines, must use "tangential none"
              //Future implementation might use this
//              if(fix_contact && ! fix_contact->handleContact(iPart,idTri,cdata.contact_history)) 
//                  continue;

              // Compute velocity of contact point on mesh
              if(vMeshC)
              {
                for(int i = 0; i < 3; i++)
                  v_wall[i] = bary[0]*vMesh[iTri][0][i] 
                            + bary[1]*vMesh[iTri][1][i] 
                            + bary[2]*vMesh[iTri][2][i];
              }
              
              cdata.i = iPart;
              cdata.deltan   = -deltan;
              cdata.delta[0] = -delta[0];
              cdata.delta[1] = -delta[1];
              cdata.delta[2] = -delta[2];
              post_force_eval_contact(cdata, mySegment,
                                      v_wall,iMesh,FixMesh_list_[iMesh],mesh,iTri);
            }
          }
        }
      }

      // clean-up contacts
      if(fix_contact) fix_contact->cleanUpContacts();
    }

}

//  **************************************************************************************
inline void FixWallGranLine::post_force_eval_contact(CollisionData & cdata,
                                                     SegmentData   & mySegment,
                                                     double * v_wall, int iMesh, 
                                                     FixMeshSurface *fix_mesh, TriMesh *mesh, int iTri)
{
  const int iPart = cdata.i;

  // but negative in distance algorithm
  cdata.r    = mySegment.radius - cdata.deltan; //corrected radius
  cdata.rsq  = cdata.r*cdata.r;
  cdata.meff = rmass_ ? rmass_[iPart] : atom->mass[atom->type[iPart]];
  cdata.area_ratio = 1.;

  cdata.computeflag = computeflag_;
  cdata.shearupdate = shearupdate_;
  cdata.jtype = atom_type_wall_;

  double force_old[3]={}, f_pw[3];

  // if force should be stored - remember old force
  if(store_force_ || stress_flag_)
    vectorCopy3D(f_[iPart],force_old);

  // add to cwl
  if(cwl_ && addflag_)
  {
      double contactPoint[3];
      double pointOnSegment[3];
      pointOnSegment[0] = mySegment.center[0] + mySegment.segmentParameter * mySegment.orientation[0];
      pointOnSegment[1] = mySegment.center[1] + mySegment.segmentParameter * mySegment.orientation[1];
      pointOnSegment[2] = mySegment.center[2] + mySegment.segmentParameter * mySegment.orientation[2];
      vectorAdd3D(pointOnSegment,cdata.delta,contactPoint);
      cwl_->add_wall_1(iMesh,mesh->id(iTri),iPart,contactPoint,v_wall);
  }

  compute_force(cdata, mySegment, v_wall); 

  // if force should be stored or evaluated
  if(store_force_ || stress_flag_)
  {
    vectorSubtract3D(f_[iPart],force_old,f_pw);

    if(store_force_)
        vectorAdd3D (wallforce_[iPart], f_pw, wallforce_[iPart]);

    if(stress_flag_ && fix_mesh->trackStress())
    {
        double delta[3];
        delta[0] = -cdata.delta[0];
        delta[1] = -cdata.delta[1];
        delta[2] = -cdata.delta[2];
        static_cast<FixMeshSurfaceStress*>(fix_mesh)->add_particle_contribution
        (
           iPart,f_pw,delta,iTri,v_wall
        );
    }
  }

  // add heat flux
  if(heattransfer_flag_)
    addHeatFlux(mesh,iPart,cdata.deltan,1.);
}

//  **************************************************************************************
void FixWallGranLine::compute_force(CollisionData & cdata, SegmentData & mySegment, double *v_wall)
{

  double **v      = atom->v;
  double **omega  = atom->omega;
  double **f      = atom->f;
  double **torque = atom->torque;
  int    *type    = atom->type;
  
  const int iPart = cdata.i;
  int itype = type[iPart];

  //0 - Compute Relative Velocity
  double uRel[3],deltaContact[3], uRelNormal;

  double distUnitVec[3];
  vectorCopy3D(cdata.delta,distUnitVec);
  vectorNormalize3D(distUnitVec); //this is the normal collision direction!
  
  for(int iDir=0;iDir<3;iDir++)
  {
    deltaContact[iDir] = mySegment.segmentParameter * mySegment.orientation[iDir] 
                       - cdata.delta[iDir]; //distance from center to contact point
  }

  vectorCross3D(omega[iPart],deltaContact,uRel);
  vectorAdd3D(uRel,v[iPart],uRel);
  vectorSubtract3D(uRel,v_wall,uRel);
  uRelNormal = -vectorDot3D(distUnitVec,uRel);
  double uRelNormalVec[3];
  vectorScalarMult3D(distUnitVec, uRelNormal, uRelNormalVec);
  vectorScalarMult3D(uRelNormalVec, -1.0);

  //1 - Pull properties
  const double kn            = kn_[itype-1][cdata.jtype-1];
  const double damping       = gamman_[itype-1][cdata.jtype-1];
  const double dampingTang   = gammat_[itype-1][cdata.jtype-1];
  const double frictionCoeff = coefficientFriction_[itype-1][cdata.jtype-1];
  const double roughness     = roughness_; 
  const double knRoughness   = kn*stiffnessFactor_;
  const double liqVisc       = fluidViscosity_;

  //2 - Compute Normal Force
  double Fn(0.0);
  // A - Region 1 (Roughness)
  if( (cdata.deltan<=0.0)  && ( -cdata.deltan <= roughness) )
       Fn = knRoughness * ( roughness + cdata.deltan );
  // B - Region 2 (Fibre)
  else if( cdata.deltan > 0.0 )
       Fn = knRoughness * roughness
          + kn          * cdata.deltan
          + damping     * uRelNormal;
  
  Fn = MathExtraLiggghts::max(0.0, Fn); //bound force to be repulsive

  // 3 - Compute Tangential Force
  double relVelTangVec[3];
  double relVelTang;
  vectorSubtract3D(uRel, uRelNormalVec , relVelTangVec);
  relVelTang = vectorMag3D(relVelTangVec);
  
  double tangUnitVec[3];
  vectorCopy3D(relVelTangVec,tangUnitVec);
  vectorNormalize3D(tangUnitVec);
  double Ft(0.0);
  Ft =  dampingTang * relVelTang; //This is a simple model, NO spring (i.e., history) interaction here!
  if (Ft > Fn*frictionCoeff)
     Ft = frictionCoeff*Fn;

  // 4 - Compute and add Normal Lubrication Force
  double Reff = mySegment.radius;
  double Flub(0.0);
  double surfSurfDist = -cdata.deltan;
  if( surfSurfDist < 2.0 * mySegment.radius ) //2.0 is a parameter here!
  {
          //TODO: check this force expression! Also: compute tangential lubrication force
          double dMaxStar = MathExtraLiggghts::max(surfSurfDist, roughness_) / Reff;
          Flub = uRelNormal * liqVisc
               * 37.699111843 //12*pi
               * Reff
               / dMaxStar;
  }
  Fn += Flub;
  
  //Apply force & torque to particle
  double normalTorque[3];
  vectorCross3D(deltaContact,distUnitVec,normalTorque);
  vectorScalarMult3D(normalTorque, Fn);

  double tangTorque[3];
  vectorCross3D(deltaContact,tangUnitVec,tangTorque);
  vectorScalarMult3D(tangTorque, Ft);

#if VERBOSEFixWallGranLine 
    #include "fix_wall_gran_line_debug.h"
#endif

  for(int iDir=0;iDir<3;iDir++)
  {
    f[iPart][iDir]      += Fn * distUnitVec[iDir]
                          +Ft * tangUnitVec[iDir];
    torque[iPart][iDir] += normalTorque[iDir] 
                          +tangTorque[iDir];
    
  }
  
  
}
