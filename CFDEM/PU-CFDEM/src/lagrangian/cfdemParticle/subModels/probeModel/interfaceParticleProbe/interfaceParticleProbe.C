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
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "interfaceParticleProbe.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"
#include "IOmanip.H"
#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interfaceParticleProbe, 0);

addToRunTimeSelectionTable
(
    probeModel,
    interfaceParticleProbe,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
interfaceParticleProbe::interfaceParticleProbe
(
    const dictionary& dict,
    cfdemCloud& sm,
    word   typeName,
    char*  logFileName
)
:
    probeModel(dict,sm,typeName,logFileName),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    twoDimensional_(false),
    depth_(1),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    pressureFieldName_(propsDict_.lookup("pressureFieldName")),
    p_(sm.mesh().lookupObject<volScalarField> (pressureFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    interpolation_(false)

   //Initialize the output stream
   {
   /* fileName probeDir = outputDirName_;
    probeDir = probeDir/U_.mesh().time().timeName();
    mkDir(probeDir);
    
    sPtrInterfaceCheck_ = new OFstream(probeDir/"InterfaceCheck");
    *sPtrInterfaceCheck_ << "#Time Velocities voidsfraction" << endl;
   */ 

    if (propsDict_.found("interpolation"))
      {
        Info << "using interpolated value of U." << endl;
        interpolation_=true;
      }

   }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interfaceParticleProbe::~interfaceParticleProbe()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interfaceParticleProbe::setForce() const
{
    label cellI=0;
    vector drag(0,0,0);
    vector Ufluid(0,0,0);
    scalar voidfraction(1);
    vector position(0,0,0);
    int numprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    OFstream* InterfaceCheckProps_;
    InterfaceCheckProps_ =  new OFstream("InterfaceCheckProps.dat");
    *InterfaceCheckProps_  << "#Time interface no-slip checking" << endl;
    
   // #include "positions.H"
  //  interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
   //Info << "positions" << positions[0] << endl;
    
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
            //for(int i=0;i<3;i++) *input_ >> positions[index][i];
            //Info << "position=" << positions << endl;
           // for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
           // {
              //  cellI = particleCloud_.cellIDs()[index][0];
                position = particleCloud_.position(index);
              //  scalar radius =  particleCloud_.radii()[index][0];
                
                //cellI = particleCloud_.cellIDs()[index][0];
                cellI = U_.mesh().findCell(position);
              //  Info << "cellI=" << cellI << endl;
                if (cellI > -1)
                {  
                   
                     if(interpolation_)
                     {
	                 //position = particleCloud_.position(index);
                        // voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                         Ufluid = UInterpolator_.interpolate(position,cellI);
                     }
                 }

      }
    
     printf("process %d \n",me);
    //Finalize and write to file
    if(Pstream::master()) //Write only if master
    {
    *InterfaceCheckProps_ << setw(IOstream::defaultPrecision() + 6) 
                << U_.mesh().time().value() << "   " 
                << Ufluid  << "   "
                << position  << "   "
                << endl;

    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

