/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Efield
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "fix.h"
#include "fix_nlocal.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "update.h"
#include "stdlib.h"
#include "mpi.h"	
#include "mpi_liggghts.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNlocal::FixNlocal(class LAMMPS *lmp, int narg, char **arg): Fix(lmp, narg, arg)
{
    int iarg = 3;	
    bool hasargs = true;

    while (iarg < narg && hasargs)
    {
	 hasargs = false;

	 if (strcmp(arg[iarg], "print_step") == 0)
	 {
		 if (iarg+2>narg) error->fix_error(FLERR, this, "not enough arguments for keyword 'print_step'");

		 this->step = atoi(arg[iarg+1]);

		 iarg += 1;
		 hasargs = true;
	 }

	 ++iarg;

    }
}

FixNlocal::~FixNlocal()
{
}


int FixNlocal::setmask()
{
    int mask = 0;
    mask |= END_OF_STEP;
    return mask;
}

void FixNlocal::end_of_step(){
	
    const int nlocal = atom->nlocal;
    const int time = update->ntimestep;
    int proc_id;
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
    
    //if( proc_id == 0 ) fprintf( screen, "FixNlocal::end_of_step() (1)  time = %d    step = %d  \n", time, this->step );
    
    if (time % this->step == 0)
    {
	    
	int nprocs = 0;
	
	MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
	MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );

	int* nproc = new int[nprocs];

	for( int i = 0; i < nprocs; ++i )
	    nproc[i] = 0;

	nproc[proc_id] = nlocal;

	MPI_Sum_Vector( nproc, nprocs, MPI_COMM_WORLD );

	if( proc_id == 0 )
	{
	     fprintf( screen, "Current timestep is %i \n", time);
	     for( int i = 0; i < nprocs; ++i )
		fprintf( screen, "Number of atoms in Processor %i: %i \n", i, nproc[i]);
	}

	delete[] nproc;
	    
    }
    
    //if( proc_id == 0 ) fprintf( screen, "FixNlocal::end_of_step() (2)  time = %d    step = %d  \n", time, step );
    
}
