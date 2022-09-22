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
#include "DriftFluxExplicit.H"
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
      defineTypeNameAndDebug(DriftFluxExplicit,0);
      addToRunTimeSelectionTable(DriftFluxModel, DriftFluxExplicit, dictionary);
    }
  }
}

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::DriftFluxExplicit::DriftFluxExplicit
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject 
)
:
   Foam::dragModels::driftFluxModels::DriftFluxModel( dict, pair, registerObject ),
   driftFluxDiffusivity
   (
       dict.subDict("driftFluxModelProps").lookup("driftFluxDiffusivity")
   ),
   filterSize_
   (
       dict.subDict("driftFluxModelProps").lookup("filterSize")   
   ),
   g_
   (
       pair.g()
   ),
   XsgsCoeff
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("XsgsCoeff"))
   ),
   YsgsCoeff
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("YsgsCoeff"))
   )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::dragModels::driftFluxModels::DriftFluxExplicit::~DriftFluxExplicit()
{}



void Foam::dragModels::driftFluxModels::DriftFluxExplicit::evolve()
{

  //fluxes
  
  dragModels::driftFluxModels::DriftFluxModel::evolve();
  Info<< "Solving for drift flux transport ... " << endl;
  
  
  volVectorField slipVelocity
  (
      pair_.continuous().U() - pair_.dispersed().U()
  );
   
  volScalarField solidVolumeFraction 
  (
        max(scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha())
  );
  
  volVectorField gradP = fvc::grad( pair_.continuous().thermo().p() );
  volVectorField gradAlpha = fvc::grad( solidVolumeFraction );
  
  
  
  volVectorField Ysgs
  (
	YsgsCoeff * solidVolumeFraction * g_
  );

  volVectorField Xsgs
  (
       XsgsCoeff * filterSize_ * filterSize_ * 
       (
          fvc::grad( gradP ) & gradAlpha
       )/max( scalar(1) - solidVolumeFraction, 1e-6 )/ pair_.continuous().rho()
  );
    
  surfaceScalarField phiUflux = fvc::interpolate( solidVolumeFraction ) * 
  				( 
				       	linearInterpolate(pair_.continuous().U()) & 												pair_.continuous().U().mesh().Sf() 
				); 
   
  //evolve drift velocity
  fvVectorMatrix driftFluxEqn
  (
      fvm::ddt( solidVolumeFraction, driftFlux_ ) 
    + fvm::div( phiUflux, driftFlux_ )
      ==
      fvm::laplacian( driftFluxDiffusivity, driftFlux_ )
    + Xsgs  
    + Ysgs
  );

  driftFluxEqn.relax();
  driftFluxEqn.solve();  
  driftFlux_.correctBoundaryConditions();
}



