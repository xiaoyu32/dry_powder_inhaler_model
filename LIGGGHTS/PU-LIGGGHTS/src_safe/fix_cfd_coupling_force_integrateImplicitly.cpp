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

   This contribution is
   Copyright 2014-     Graz University of Technology, IPPT (radl@tugraz.at)
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_force_integrateImplicitly.h"
#include "fix_property_particle.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

enum{NONE,DIPOLE};
/* ---------------------------------------------------------------------- */

FixCfdCouplingForceIntegrateImplicitly::FixCfdCouplingForceIntegrateImplicitly(LAMMPS *lmp, int narg, char **arg) :
    FixCfdCouplingForce(lmp,narg,arg),
    CNalpha_(0.5),
    CAddRhoFluid_(0.0),
    onePlusCAddRhoFluid_(1.0),
    fix_Ksl_(0),
    fix_uf_(0)
{
    int iarg = 3;

    // process extra keywords
    extra = NONE;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"CrankNicolson") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'CrankNicholson'");
            iarg++;
            CNalpha_ = atof(arg[iarg]);
            if(CNalpha_<0 || CNalpha_>1)
                error->fix_error(FLERR,this,"incorrect choice for 'CrankNicholson': setting CNalpha_<0 or CNalpha_>1 is not appropriate");

            fprintf(screen,"cfd_coupling_foce_implicit will use Crank-Nicholson scheme with %f\n", CNalpha_);
            iarg++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"CAddRhoFluid") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'CAddRhoFluid'");
            iarg++;
            CAddRhoFluid_        = atof(arg[iarg]);
            onePlusCAddRhoFluid_ = 1.0 + CAddRhoFluid_;
            fprintf(screen,"cfd_coupling_foce_implicit will consider added mass with CAddRhoFluid = %f\n",
                    CAddRhoFluid_);
            iarg++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"update") == 0)
        {
          if (iarg+2 > narg) error->all(FLERR,"Illegal fix nve/sphere command");
          if (strcmp(arg[iarg+1],"dipole") == 0) extra = DIPOLE;
          else error->all(FLERR,"Illegal fix nve/sphere command");
          iarg += 2;
        }
    }
  nevery = 1;
}

/* ---------------------------------------------------------------------- */
FixCfdCouplingForceIntegrateImplicitly::~FixCfdCouplingForceIntegrateImplicitly(){}

/* ---------------------------------------------------------------------- */
int FixCfdCouplingForceIntegrateImplicitly::setmask()
{
  int mask = 0;

  //A - NVE Relevant
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;

  //B - COUPLING Relevant
  mask |= POST_FORCE;
  mask |= END_OF_STEP;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceIntegrateImplicitly::post_create()
{
    // do mother class init w/o dragforce
    FixCfdCouplingForce::post_create();

    // register Ksl
    if(!fix_Ksl_)
    {
        const char* fixarg[9];
        fixarg[0]="Ksl";
        fixarg[1]="all";
        fixarg[2]="property/particle";
        fixarg[3]="Ksl";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_Ksl_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register uf
    if(!fix_uf_)
    {
        const char* fixarg[11];
        fixarg[0]="uf";
        fixarg[1]="all";
        fixarg[2]="property/particle";
        fixarg[3]="uf";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_uf_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceIntegrateImplicitly::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_Ksl_) modify->delete_fix("Ksl");
    if(unfixflag && fix_uf_) modify->delete_fix("uf");
}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingForceIntegrateImplicitly::init()
{
  //A - Init NVE part
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;

  // check that all particles are finite-size spheres
  // no point particles allowed
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nve/sphere requires extended particles");

  //B - Init coupling part
  FixCfdCouplingForce::init();
  fix_coupling_->add_pull_property("Ksl","scalar-atom");
  fix_coupling_->add_pull_property("uf","vector-atom");
}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingForceIntegrateImplicitly::initial_integrate(int vflag)
{
  double dtirotate,msq,scale;
  double g[3];

  double **x       = atom->x;
  double **v       = atom->v;
  double **f       = atom->f;
  double **omega   = atom->omega;
  double **torque  = atom->torque;
  double *radius   = atom->radius;
  double *rmass    = atom->rmass;
  int    *mask     = atom->mask;

  double *Ksl = fix_Ksl_->vector_atom;
  double **uf = fix_uf_->array_atom;

  int nlocal = atom->nlocal;

  double frc[3];

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate;
  if (domain->dimension == 2) dtfrotate = dtf / 0.5; // for discs the formula is I=0.5*Mass*Radius^2
  else dtfrotate  = dtf / INERTIA;

  // update 1/2 step for v and omega, and full step for  x for all particles
  // d_omega/dt = torque / inertia
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

/*      printf("initial_integrate: f: %g %g %g, v: %g %g %g, Ksl: %g, uf: %g %g %g \n",
             f[i][0],f[i][1],f[i][2],
             v[i][0],v[i][1],v[i][2],
             Ksl[i],
             uf[i][0],uf[i][1],uf[i][2]
             );*/

      //Implicit 1/2 step of implicit velocity update, frc is also half force!
      implicitVelocityUpdate( dtf, rmass[i],
                              v[i], f[i], Ksl[i], uf[i],
                              frc
                            );

      // add up forces for post-proc
      vectorAdd3D(dragforce_total,frc,dragforce_total);

      x[i][0] += dtv  * v[i][0];
      x[i][1] += dtv  * v[i][1];
      x[i][2] += dtv  * v[i][2];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
  }

  // update mu for dipoles
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length
  if (extra == DIPOLE) {
    double **mu = atom->mu;
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
  }
}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingForceIntegrateImplicitly::post_force(int)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **dragforce = fix_dragforce_->array_atom;

  vectorZeroize3D(dragforce_total);

  // add dragforce to force vector
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        // add NON-Implicit fluid-particle forces
        vectorAdd3D(f[i],dragforce[i],f[i]);

        // add up forces for post-proc
        vectorAdd3D(dragforce_total,dragforce[i],dragforce_total);
    }
  }

}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingForceIntegrateImplicitly::final_integrate()
{
  double dtirotate;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int    *mask = atom->mask;
  int     nlocal = atom->nlocal;

  double *Ksl = fix_Ksl_->vector_atom;
  double **uf = fix_uf_->array_atom;

  double frc[3];

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  double dtfrotate;
  if (domain->dimension == 2) dtfrotate = dtf / 0.5; // for discs the formula is I=0.5*Mass*Radius^2
  else dtfrotate  = dtf / INERTIA;

  // update 1/2 stef for v and omega for all particles
  // d_omega/dt = torque / inertia
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

/*      printf("final_integrate: f: %g %g %g, v: %g %g %g, Ksl: %g, uf: %g %g %g \n",
             f[i][0],f[i][1],f[i][2],
             v[i][0],v[i][1],v[i][2],
             Ksl[i],
             uf[i][0],uf[i][1],uf[i][2]
             );*/

      //Implicit 1/2 step of implicit velocity update
      implicitVelocityUpdate( dtf, rmass[i],
                              v[i], f[i], Ksl[i], uf[i],
                              frc
                            );

      // add up forces for post-proc
      vectorAdd3D(dragforce_total,frc,dragforce_total);

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingForceIntegrateImplicitly::end_of_step()
{
}

/* --------------------------------------------------------------------- */
void FixCfdCouplingForceIntegrateImplicitly::implicitVelocityUpdate
            (
                double dtf, double rmass,
                double *v, double *f, double Ksl, double *uf,
                double *frc
            )
{
      double vN32[3], deltaU, dtfm, KslMDeltaT;

      dtfm       = dtf / ( rmass*onePlusCAddRhoFluid_ );
      KslMDeltaT = Ksl * dtfm;

      for(int dirI=0;dirI<3;dirI++)
      {
            //calculate new velocity
            vN32[dirI] = (  v[dirI]
                          + f[dirI] * dtfm
                          + KslMDeltaT
                           * (   uf[dirI]
                              - (1.0-CNalpha_)*v[dirI]
                             )
                         )
                         /
                         (1.0+KslMDeltaT*CNalpha_);

            //calculate velocity difference sand force
            deltaU    =  uf[dirI]
                      - (
                            (1.0-CNalpha_)*v[dirI]
                           +     CNalpha_ *vN32[dirI]
                        );
            frc[dirI]  = Ksl * deltaU * 0.5;  //1//2 of IMPLICIT drag force

            //update the particle velocity
            v[dirI] = vN32[dirI];  //update velocity for a half step!
      }
}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingForceIntegrateImplicitly::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingForceIntegrateImplicitly::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */
void FixCfdCouplingForceIntegrateImplicitly::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
