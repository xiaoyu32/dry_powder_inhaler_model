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
#include "OzelIJMF.H"
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
      defineTypeNameAndDebug(OzelIJMF,0);
      addToRunTimeSelectionTable(DriftFluxModel, OzelIJMF, dictionary);
    }
  }
}

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::OzelIJMF::OzelIJMF
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject 
)
:
   Foam::dragModels::driftFluxModels::DriftFluxModel( dict, pair, registerObject )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::dragModels::driftFluxModels::OzelIJMF::~OzelIJMF()
{}

void Foam::dragModels::driftFluxModels::OzelIJMF::evolve()
{
  dragModels::driftFluxModels::DriftFluxModel::evolve();
  Info<< "Ozel-IJFM drift flux model... " << endl;
  
  volVectorField slipVelocity
  (
      pair_.continuous().U() - pair_.dispersed().U()
  );
   
  volScalarField solidVolumeFraction 
  (
        max(scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha())
  );

  scalar alphaMax(0.64);
  scalar Kbb(3.33);
  scalar func_f(1);
  scalar constH1(0.2);
  scalar constH2(-1.88);
  scalar constH3(5.16);
	
  volScalarField func_h = -  tanh(solidVolumeFraction/constH1)
                            *sqrt(solidVolumeFraction/alphaMax)
                            *sqr(1.-solidVolumeFraction/alphaMax)
                            *(1.+constH2*solidVolumeFraction/alphaMax+constH3*sqr(solidVolumeFraction/alphaMax));
 
  //evolve drift velocity
  driftFlux_ = Kbb*func_f*func_h*solidVolumeFraction*slipVelocity;
  driftFlux_.correctBoundaryConditions();
}

