/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "DriftFluxModel.H"
#include "ConstantDriftFluxCorr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
//#include "fvOption.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace dragModels
  {
    namespace driftFluxModels
    {
      defineTypeNameAndDebug(ConstantDriftFluxCorr,0);
      addToRunTimeSelectionTable(DriftFluxModel, ConstantDriftFluxCorr, dictionary);
    }
  }
}

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::ConstantDriftFluxCorr::ConstantDriftFluxCorr
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject 
)
:
   Foam::dragModels::driftFluxModels::DriftFluxModel( dict, pair, registerObject ),
   constantCorr_ 
   (
       readScalar(dict.subDict("ConstantDriftFluxCorrProps").lookup("correctionConstant"))
   )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::dragModels::driftFluxModels::ConstantDriftFluxCorr::~ConstantDriftFluxCorr()
{}

void Foam::dragModels::driftFluxModels::ConstantDriftFluxCorr::evolve()
{
  dragModels::driftFluxModels::DriftFluxModel::evolve();
  Info<< "Constant drift flux model... " << endl;
  
  volVectorField slipVelocity
  (
      pair_.continuous().U() - pair_.dispersed().U()
  );
   
  volScalarField solidVolumeFraction 
  (
        max(scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha())
  );
  
  //evolve drift velocity
  driftFlux_ = constantCorr_*solidVolumeFraction*slipVelocity;
  driftFlux_.correctBoundaryConditions();
}

