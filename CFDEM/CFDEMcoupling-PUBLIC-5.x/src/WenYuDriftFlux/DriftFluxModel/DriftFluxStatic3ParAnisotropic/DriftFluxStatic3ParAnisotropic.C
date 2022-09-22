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
#include "DriftFluxStatic3ParAnisotropic.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
//#include "fvOption.H"
#include "phasePair.H"
#include "extendedCentredCellToCellStencil.H"
#include "MeshObject.H"
#include "keras_model.H"
#include <vector>



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace dragModels
  {
    namespace driftFluxModels
    {
      defineTypeNameAndDebug(DriftFluxStatic3ParAnisotropic,0);
      addToRunTimeSelectionTable(DriftFluxModel, DriftFluxStatic3ParAnisotropic, dictionary);
    }
  }
}

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::DriftFluxStatic3ParAnisotropic::DriftFluxStatic3ParAnisotropic
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject 
)
:
   Foam::dragModels::driftFluxModels::DriftFluxModel( dict, pair, registerObject ),
   mesh_( pair.dispersed().mesh() ),
   gradP
   (
	IOobject
	(
            "gradP",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh(),
	dimensionedVector("gradP",dimensionSet(1,-2,-2,0,0),Foam::vector(0,0,0))
   ),  
   gradAlphaP
   (
	IOobject
	(
            "gradAlphaP",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh(),
	dimensionedVector("gradAlphaP",dimensionSet(0,-1,0,0,0),Foam::vector(0,0,0))
   ),  
   DFXnnModel_
   (
       "DFXkerasParameters.nnet",
       false
   ),
   DFYnnModel_
   (
       "DFYkerasParameters.nnet",
       false
   ),
   DFZnnModel_
   (
       "DFZkerasParameters.nnet",
       false
   ),
   uTerminal
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("terminalVelocity"))
   ),
   lStar
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("characteristicLength"))
   ),
   alphaMax
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("MaximumSolidVolumeFraction"))
   ),
   filterSize_
   (
       dict.subDict("driftFluxModelProps").lookup("filterSize")   
   ),
   g_
   (
       pair.g()
   ), 
   uSlip
   (
      pair_.continuous().U() - pair_.dispersed().U()
   ),
   uSlipVdrift
   ( 
      pair_.continuous().U() - pair_.dispersed().U() + driftFlux_
   ),
   solidVolumeFraction 
   (
       max( scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha() )
   ),
   rhoParticle_
   (
	pair_.dispersed().rho()[0]
    ),
   rhoFluid_
   (
	pair_.continuous().rho()[0]
    ),
   DFNFeatures
   (
       int( readScalar(dict.subDict("driftFluxModelProps").lookup("DF_NumberOfFeatures")) )  
   ),
   DFXmeans_(NULL),
   DFXvars_(NULL),
   DFYmeans_(NULL),
   DFYvars_(NULL),
   DFZmeans_(NULL),
   DFZvars_(NULL)
{
   
  //Info<<"Number of features: "<< DFNFeatures <<endl;
   
  DFXmeans_ = new double[DFNFeatures];
  DFXvars_ = new double[DFNFeatures];
  DFYmeans_ = new double[DFNFeatures];
  DFYvars_ = new double[DFNFeatures];
  DFZmeans_ = new double[DFNFeatures];
  DFZvars_ = new double[DFNFeatures];
  
     
    std::ifstream reader;

    reader.open( "DFX_mean.csv" );
    
    if( !reader.is_open() )
    {
        Info<<"Failed to open file for DF_mean!"<<endl;
	
	for( int i = 0; i < DFNFeatures; ++i )
	    DFXmeans_[i] = 0.0;
	
    }else
    {
	for( int i = 0; i < DFNFeatures; ++i )
            reader>>DFXmeans_[i];
    }
    
    reader.close();

    reader.open( "DFX_std.csv" );

    if( !reader.is_open() )
    {
        Info<<"Failed to open file for DF_std!"<<endl;
	
	for( int i = 0; i < DFNFeatures; ++i )
	    DFXvars_[i] = 1.0;
	
    }else{
    
	for( int i = 0; i < DFNFeatures; ++i )
            reader>>DFXvars_[i];
	    
    }
    
    reader.close();
    reader.open( "DFY_mean.csv" );
    
    if( !reader.is_open() )
    {
        Info<<"Failed to open file for DF_mean!"<<endl;
	
	for( int i = 0; i < DFNFeatures; ++i )
	    DFYmeans_[i] = 0.0;
	
    }else
    {
	for( int i = 0; i < DFNFeatures; ++i )
            reader>>DFYmeans_[i];
    }
    
    reader.close();

    reader.open( "DFY_std.csv" );

    if( !reader.is_open() )
    {
        Info<<"Failed to open file for DF_std!"<<endl;
	
	for( int i = 0; i < DFNFeatures; ++i )
	    DFYvars_[i] = 1.0;
	
    }else{
    
	for( int i = 0; i < DFNFeatures; ++i )
            reader>>DFYvars_[i];
	    
    }
    
    reader.close();
    reader.open( "DFZ_mean.csv" );
    
    if( !reader.is_open() )
    {
        Info<<"Failed to open file for DF_mean!"<<endl;
	
	for( int i = 0; i < DFNFeatures; ++i )
	    DFZmeans_[i] = 0.0;
	
    }else
    {
	for( int i = 0; i < DFNFeatures; ++i )
            reader>>DFZmeans_[i];
    }
    
    reader.close();

    reader.open( "DFZ_std.csv" );

    if( !reader.is_open() )
    {
        Info<<"Failed to open file for DF_std!"<<endl;
	
	for( int i = 0; i < DFNFeatures; ++i )
	    DFZvars_[i] = 1.0;
	
    }else{
    
	for( int i = 0; i < DFNFeatures; ++i )
            reader>>DFZvars_[i];
	    
    }
    
    reader.close();

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::dragModels::driftFluxModels::DriftFluxStatic3ParAnisotropic::~DriftFluxStatic3ParAnisotropic()
{
  if( DFXmeans_ )
  {
     delete[] DFXmeans_;
  }
  if( DFXvars_ )
  {
     delete[] DFXvars_;
  }
  if( DFYmeans_ )
  {
     delete[] DFYmeans_;
  }
  if( DFYvars_ )
  {
     delete[] DFYvars_;
  }
  if( DFZmeans_ )
  {
     delete[] DFZmeans_;
  }
  if( DFZvars_ )
  {
     delete[] DFZvars_;
  }
}



void Foam::dragModels::driftFluxModels::DriftFluxStatic3ParAnisotropic::evolve()
{

  //fluxes
  
  dragModels::driftFluxModels::DriftFluxModel::evolve();
  Info<< "Solving for drift flux transport ... " << endl;
  solidVolumeFraction = max( scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha() );
  uSlip = pair_.continuous().U() - pair_.dispersed().U();
  //  gradAlphaP = fvc::grad( solidVolumeFraction );
  gradP = fvc::grad( pair_.continuous().thermo().p());
  //driftFlux_min = -solidVolumeFraction*uSlip;
  //  tmp<volVectorField> gradGradPz = fvc::grad(gradP.component(2));
  //tmp<volScalarField> laplaceGradPz = fvc::laplacian(gradP.component(2));
  //tmp<volScalarField> laplaceAlphaP = fvc::laplacian(solidVolumeFraction);

  /******* Neural Network for DriftFlux********/


  //const extendedCentredCellToCellStencil& stencil = this->stencil();
  
  #include "DF.H"



  
  forAll( driftFlux_, iCell )
  {    
     if( solidVolumeFraction[iCell] <= 1e-3 )
     {
         driftFlux_[iCell] = Foam::vector( 0.0, 0.0, 0.0 );
     }
  }
  driftFlux_.correctBoundaryConditions();

}

