/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

    two way DEM-CFD coupling via MPI

Class
    twoWayMPIpu

SourceFiles
    twoWayMPIpu.C

\*---------------------------------------------------------------------------*/

#ifndef twoWayMPIpu_H
#define twoWayMPIpu_H

#include "dataExchangeModel.H"
#include "liggghtsCommandModel.H"
#include "OFstream.H"
#include "sys/stat.h"
#include "pair.h"
#include "force.h"
#include "forceModel.H"

//=================================//
//LAMMPS/LIGGGHTS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include <lammps.h>         // these are LAMMPS include files
#include <input.h>
#include <atom.h>
#include <library.h>
#include <error.h>
#include <library_cfd_coupling.h>
#include <update.h>
//=================================//

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/


class twoWayMPIpu
:
    public dataExchangeModel
{
private:

  // private data
    dictionary propsDict_;

    // variables
    int me;

    int nprocs;

    int liggghts;

    MPI_Comm comm_liggghts;

    LAMMPS_NS::LAMMPS *lmp;

  // private member functions

    //char* wordToChar(word&) const; // mostafa change this to const
    char* wordToChar(word&) const;

    mutable int bufferSizeForAlloc_;
    
    const static double asciiPrecission = 1e-6;
    
public:

    //- Runtime type information
    TypeName("twoWayMPIpu");


    // Constructors

        //- Construct from components
        twoWayMPIpu
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~twoWayMPIpu();


    // Member Functions
        void getData
        (
            word name,
            word type,
            double ** const& field,
            label step
        ) const;

        void getData
        (
            word name,
            word type,
            int ** const& field,
            label step
        ) const;

        void giveData
        (
            word name,
            word type,
            double ** const& field,
            const char* datatype
        ) const;
	
	virtual void* getLmp() const;
	
        //============
        // double **
        void allocateArray(double**&, double, int, int) const;
        void allocateArray(double**&, double, int,const char* ="nparticles") const;
        void destroy(double**,int) const;
        //============
        // int **
        void allocateArray(int**&, int, int, int) const;
        void allocateArray(int**&, int, int,const char* ="nparticles") const;
        void destroy(int**,int) const;
        //==============
        // double *
        void destroy(double*) const;
        //============
        // int **
        void destroy(int*) const;
        //==============

	void write_boundary_info( cfdemCloud& );
	
        bool couple() const;

        int getNumberOfParticles() const;
        int getNumberOfClumps() const;

        word myType() const{return typeName; };

//_AO 06/02/2014
//        void setCG() const { particleCloud_.forceM(0).setCG(lmp->force->cg()); };
	void setCG() const { particleCloud_.setCG(lmp->force->cg()); };
//_AO 06/02/2014

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
