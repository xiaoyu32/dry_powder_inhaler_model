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
#include <cstdio>
#include <cstring>
#include "fix_property_global.h"
#include "fix_nve_sphere_parcel.h"
#include "fix_nve_sphere.h"
#include "atom.h"
#include "domain.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "math_vector.h"
#include "math_extra.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

enum{NONE,DIPOLE};
enum{NODLM,DLM};

/* ---------------------------------------------------------------------- */

FixNVESphereParcel::FixNVESphereParcel(LAMMPS *lmp, int narg, char **arg) :
  FixNVESphere(lmp, narg, arg),
  fix_nparcel( NULL ),
  nparcel( NULL )
{}

FixNVESphereParcel::~FixNVESphereParcel()
{
   if( nparcel ) delete[] nparcel;
}

/* ---------------------------------------------------------------------- */

void FixNVESphereParcel::init()
{
  FixNVESphere::init();
  int max_type = atom->get_properties()->max_type();
  fix_nparcel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("nparcel","property/particle","peratomtype",max_type,0,style));
  
  if( !fix_nparcel ) error->all(FLERR,"Integration style 'nve/sphere/parcel' requires property/particle fix 'nparcel'!");
  
  nparcel = new double[max_type];
  
  for( int i = 0; i < max_type; ++i )
  {
      nparcel[i] = fix_nparcel->compute_vector(i);
      
      if( nparcel[i] <= 0 ) error->all(FLERR,"Number of particles in a parcel has to be a positive number!");
         
  }    
      
  
}

/* ---------------------------------------------------------------------- */

void FixNVESphereParcel::initial_integrate(int /*vflag*/)
{
  double dtfm,dtirotate,msq,scale,s2,inv_len_mu;
  double g[3];
  vector w, w_temp, a;
  matrix Q, Q_temp, R;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  
  int *type = atom->type;
  
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / inertia;

  // update v,x,omega for all particles
  // d_omega/dt = torque / inertia
    
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
    
      //std::cout<<"DEM type: "<<atom->type[i]<<std::endl;
      //fprintf( screen, "DEM type = %d \n", atom->type[i] );
       
      dtfm = dtf / ( rmass[i] * nparcel[type[i]-1]  );
      
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i] * nparcel[type[i]-1] );
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
  }

  // update mu for dipoles

  if (extra == DIPOLE) {
    double **mu = atom->mu;
    if (dlm == NODLM) {

      // d_mu/dt = omega cross mu
      // renormalize mu to dipole length

      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          if (mu[i][3] > 0.0) {
            g[0] = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
            g[1] = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
            g[2] = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
            msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
            scale = mu[i][3]/sqrt(msq);
            mu[i][0] = g[0]*scale;
            mu[i][1] = g[1]*scale;
            mu[i][2] = g[2]*scale;
          }
    } else {

      // integrate orientation following Dullweber-Leimkuhler-Maclachlan scheme

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit && mu[i][3] > 0.0) {

          // Construct Q from dipole:
          // Q is the rotation matrix from space frame to body frame
          // i.e. v_b = Q.v_s

          // define mu to lie along the z axis in the body frame
          // take the unit dipole to avoid getting a scaling matrix

          inv_len_mu = 1.0/mu[i][3];
          a[0] = mu[i][0]*inv_len_mu;
          a[1] = mu[i][1]*inv_len_mu;
          a[2] = mu[i][2]*inv_len_mu;

          // v = a x [0 0 1] - cross product of mu in space and body frames
          // s = |v|
          // c = a.[0 0 1] = a[2]
          // vx = [ 0    -v[2]  v[1]
          //        v[2]  0    -v[0]
          //       -v[1]  v[0]  0    ]
          // then
          // Q = I + vx + vx^2 * (1-c)/s^2

          s2 = a[0]*a[0] + a[1]*a[1];
          if (s2 != 0.0){ // i.e. the vectors are not parallel
            scale = (1.0 - a[2])/s2;

            Q[0][0] = 1.0 - scale*a[0]*a[0];
            Q[0][1] = -scale*a[0]*a[1];
            Q[0][2] = -a[0];
            Q[1][0] = -scale*a[0]*a[1];
            Q[1][1] = 1.0 - scale*a[1]*a[1];
            Q[1][2] = -a[1];
            Q[2][0] = a[0];
            Q[2][1] = a[1];
            Q[2][2] = 1.0 - scale*(a[0]*a[0] + a[1]*a[1]);
          } else { // if parallel then we just have I or -I
            Q[0][0] = 1.0/a[2];  Q[0][1] = 0.0;       Q[0][2] = 0.0;
            Q[1][0] = 0.0;       Q[1][1] = 1.0/a[2];  Q[1][2] = 0.0;
            Q[2][0] = 0.0;       Q[2][1] = 0.0;       Q[2][2] = 1.0/a[2];
          }

          // Local copy of this particle's angular velocity (in space frame)
          w[0] = omega[i][0]; w[1] = omega[i][1]; w[2] = omega[i][2];

          // Transform omega into body frame: w_temp= Q.w
          matvec(Q,w,w_temp);

          // Construct rotation R1
          BuildRxMatrix(R, dtf/force->ftm2v*w_temp[0]);

          // Apply R1 to w: w = R.w_temp
          matvec(R,w_temp,w);

          // Apply R1 to Q: Q_temp = R^T.Q
          transpose_times3(R,Q,Q_temp);

          // Construct rotation R2
          BuildRyMatrix(R, dtf/force->ftm2v*w[1]);

          // Apply R2 to w: w_temp = R.w
          matvec(R,w,w_temp);

          // Apply R2 to Q: Q = R^T.Q_temp
          transpose_times3(R,Q_temp,Q);

          // Construct rotation R3
          BuildRzMatrix(R, 2.0*dtf/force->ftm2v*w_temp[2]);

          // Apply R3 to w: w = R.w_temp
          matvec(R,w_temp,w);

          // Apply R3 to Q: Q_temp = R^T.Q
          transpose_times3(R,Q,Q_temp);

          // Construct rotation R4
          BuildRyMatrix(R, dtf/force->ftm2v*w[1]);

          // Apply R4 to w: w_temp = R.w
          matvec(R,w,w_temp);

          // Apply R4 to Q: Q = R^T.Q_temp
          transpose_times3(R,Q_temp,Q);

          // Construct rotation R5
          BuildRxMatrix(R, dtf/force->ftm2v*w_temp[0]);

          // Apply R5 to w: w = R.w_temp
          matvec(R,w_temp,w);

          // Apply R5 to Q: Q_temp = R^T.Q
          transpose_times3(R,Q,Q_temp);

          // Transform w back into space frame w_temp = Q^T.w
          transpose_matvec(Q_temp,w,w_temp);
          omega[i][0] = w_temp[0];
          omega[i][1] = w_temp[1];
          omega[i][2] = w_temp[2];

          // Set dipole according to updated Q: mu = Q^T.[0 0 1] * |mu|
          mu[i][0] = Q_temp[2][0] * mu[i][3];
          mu[i][1] = Q_temp[2][1] * mu[i][3];
          mu[i][2] = Q_temp[2][2] * mu[i][3];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereParcel::final_integrate()
{
  double dtfm,dtirotate;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / inertia;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia
   
  double rke = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / ( rmass[i] * nparcel[type[i]-1] );
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]* nparcel[type[i]-1]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
      rke += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] +
              omega[i][2]*omega[i][2])*radius[i]*radius[i]*rmass[i];
    }

}
