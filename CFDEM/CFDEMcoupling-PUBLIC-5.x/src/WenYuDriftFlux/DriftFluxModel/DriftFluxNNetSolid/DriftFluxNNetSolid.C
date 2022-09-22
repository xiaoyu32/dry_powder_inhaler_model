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
#include "DriftFluxNNetSolid.H"
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
      defineTypeNameAndDebug(DriftFluxNNetSolid,0);
      addToRunTimeSelectionTable(DriftFluxModel, DriftFluxNNetSolid, dictionary);
    }
  }
}

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::DriftFluxNNetSolid::DriftFluxNNetSolid
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject 
)
:
   Foam::dragModels::driftFluxModels::DriftFluxModel( dict, pair, registerObject ),
   mesh_( pair.dispersed().mesh() ),
   Xsgs
   (
	IOobject
	(
            "Xsgs",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh(),
	dimensionedVector("Xsgs",dimensionSet(0,1,-2,0,0),Foam::vector(0,0,0))
   ),
   Ysgs
   (
	IOobject
	(
            "Ysgs",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh(),
	dimensionedVector("Ysgs",dimensionSet(0,1,-2,0,0),Foam::vector(0,0,0))
   ),  
   sumXY
   (
	IOobject
	(
            "sumXY",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh(),
	dimensionedVector("sumXY",dimensionSet(0,1,-2,0,0),Foam::vector(0,0,0))
   ), 
   driftFlux_min
   (
    (max( scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha() ))*(pair_.continuous().U() - pair_.dispersed().U())
    ),

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
   XnnModel_
   (
       "XkerasParameters.nnet",
       false
   ),
   YnnModel_
   (
       "YkerasParameters.nnet",
       false
   ),
   XBnnModel_
   (
       "XBkerasParameters.nnet",
       false
   ),
   driftFluxDiffusivity
   (
       dict.subDict("driftFluxModelProps").lookup("driftFluxDiffusivity")
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
   /*
   gradP
   (
       fvc::grad( pair_.continuous().thermo().p())
   ),
   gradAlphaP 
   ( 
       fvc::grad( solidVolumeFraction )
   ), 
   */
   XNFeatures
   (
       int( readScalar(dict.subDict("driftFluxModelProps").lookup("X_NumberOfFeatures")) )  
   ),
   XBNFeatures
   (
       int( readScalar(dict.subDict("driftFluxModelProps").lookup("XB_NumberOfFeatures")) )  
   ),
   YNFeatures
   (
       int( readScalar(dict.subDict("driftFluxModelProps").lookup("Y_NumberOfFeatures")) )
   ),
   Ymin
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("Ymin"))
    ),
   Ymax
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("Ymax"))
    ) ,
   Xmin
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("Xmin"))
    ),
   Xmax
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("Xmax"))
    ),
   Vdmax
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("Vdmax"))
    ),
   Vdmin
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("Vdmin"))
    ),
   sumXYmax
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("sumXYmax"))
    ),
   sumXYmin
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("sumXYmin"))
    ),
    Xmeans_(NULL),
    Xvars_(NULL),
    XBmeans_(NULL),
    XBvars_(NULL),
    Ymeans_(NULL),
    Yvars_(NULL)
  
{
   
  Info<<"Number of features: "<< XNFeatures <<" / "<<YNFeatures<<endl;
   
  Xmeans_ = new double[XNFeatures];
  Xvars_ = new double[XNFeatures];
  Ymeans_ = new double[YNFeatures];
  Yvars_ = new double[YNFeatures];
  XBmeans_ = new double[XBNFeatures];
  XBvars_ = new double[XBNFeatures];
  
     
    std::ifstream reader;

    reader.open( "X_mean.csv" );
    
    if( !reader.is_open() )
    {
        Info<<"Failed to open file for X_mean!"<<endl;
	
	for( int i = 0; i < XNFeatures; ++i )
	    Xmeans_[i] = 0.0;
	
    }else{
	for( int i = 0; i < XNFeatures; ++i )
            reader>>Xmeans_[i];
    }
    
    reader.close();
    reader.open( "XB_mean.csv" );
    
    if( !reader.is_open() )
    {
        Info<<"Failed to open file for XB_mean!"<<endl;
	
	for( int i = 0; i < XBNFeatures; ++i )
	    XBmeans_[i] = 0.0;
	
    }else{
	for( int i = 0; i < XBNFeatures; ++i )
            reader>>XBmeans_[i];
    }
    
    reader.close();

    reader.open( "X_std.csv" );

    if( !reader.is_open() )
    {
        Info<<"Failed to open file for X_std!"<<endl;
	
	for( int i = 0; i < XNFeatures; ++i )
	    Xvars_[i] = 1.0;
	
    }else{
    
	for( int i = 0; i < XNFeatures; ++i )
            reader>>Xvars_[i];
	    
    }
    
    reader.close();

    reader.open( "XB_std.csv" );

    if( !reader.is_open() )
    {
        Info<<"Failed to open file for XB_std!"<<endl;
	
	for( int i = 0; i < XBNFeatures; ++i )
	    XBvars_[i] = 1.0;
	
    }else{
    
	for( int i = 0; i < XBNFeatures; ++i )
            reader>>XBvars_[i];
	    
    }
    
    reader.close();



    reader.open( "Y_mean.csv" );
    
    if( !reader.is_open() )
    {
        Info<<"Failed to open file for Y_mean!"<<endl;
	
	for( int i = 0; i < YNFeatures; ++i )
	    Ymeans_[i] = 0.0;
	
    }else{
	for( int i = 0; i < YNFeatures; ++i )
            reader>>Ymeans_[i];
    }
    
    reader.close();



    reader.open( "Y_std.csv" );

    if( !reader.is_open() )
    {
        Info<<"Failed to open file for Y_std!"<<endl;
	
	for( int i = 0; i < YNFeatures; ++i )
	    Yvars_[i] = 1.0;
	
    }else{
    
	for( int i = 0; i < YNFeatures; ++i )
            reader>>Yvars_[i];
	    
    }
    
    reader.close();
      
      

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::dragModels::driftFluxModels::DriftFluxNNetSolid::~DriftFluxNNetSolid()
{
  if( Xmeans_ )
  {
     delete[] Xmeans_;
  }
  if( Xvars_ )
  {
     delete[] Xvars_;
  }
  if( XBmeans_ )
  {
     delete[] XBmeans_;
  }
  if( XBvars_ )
  {
     delete[] XBvars_;
  }    
  if( Ymeans_ )
  {
     delete[] Ymeans_;
  }
    if( Yvars_ )
  {
     delete[] Yvars_;
  }

}



void Foam::dragModels::driftFluxModels::DriftFluxNNetSolid::evolve()
{

  //fluxes
  
  dragModels::driftFluxModels::DriftFluxModel::evolve();
  Info<< "Solving for drift flux transport ... " << endl;
  solidVolumeFraction = max( scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha() );
  uSlip = pair_.continuous().U() - pair_.dispersed().U();
  gradAlphaP = fvc::grad( solidVolumeFraction );
  gradP = fvc::grad( pair_.continuous().thermo().p());
  driftFlux_min = -solidVolumeFraction*uSlip;

  /******* Neural Network for Ysgs********/


  const extendedCentredCellToCellStencil& stencil = this->stencil();
  
  #include "Ysgs_z.H"

  /******** Neural Network for Xsgs *************/

  
  #include "Xsgs_z.H"

  /*
  surfaceScalarField phiUflux = fvc::interpolate( solidVolumeFraction ) * 
  				( 
				       	linearInterpolate(pair_.continuous().U()) & 									pair_.continuous().U().mesh().Sf() 
				); 
  */
  surfaceScalarField phiUflux =  linearInterpolate(pair_.dispersed().U()) & 									pair_.dispersed().U().mesh().Sf(); 
  //evolve drift velocity
  

//with diffusion
  fvVectorMatrix driftFluxEqn
  (
      fvm::ddt(  driftFlux_ ) 
    + fvm::div( phiUflux, driftFlux_ )
      ==
      fvm::laplacian( driftFluxDiffusivity, driftFlux_ )
      //+ sumXY
      + Xsgs/(1.0-solidVolumeFraction)  
      + Ysgs
  );
  


  /*
  //without diffusion
  fvVectorMatrix driftFluxEqn
  (
   //fvm::ddt( solidVolumeFraction, driftFlux_ ) 
   fvm::ddt( driftFlux_ )
    + fvm::div( phiUflux, driftFlux_ )
      ==
        Xsgs/(1.0-solidVolumeFraction)  
      + Ysgs
  );
  */

  driftFluxEqn.relax();
  driftFluxEqn.solve();  
  
  forAll( driftFlux_, iCell )
  {    
     if( solidVolumeFraction[iCell] <= 1e-3 )
     {
         driftFlux_[iCell] = Foam::vector( 0.0, 0.0, 0.0 );
     }
     
     //driftFlux_[iCell].z() = max(min(driftFlux_[iCell].z(),Vdmax),Vdmin);
     else
     {
       
	   if (uSlip[iCell].z()>=0)
	     {
	       if (driftFlux_[iCell].z() < (0.93*driftFlux_min[iCell].z())) driftFlux_[iCell].z() = 0.93*driftFlux_min[iCell].z();
	       if (driftFlux_[iCell].z() > (-0.05*driftFlux_min[iCell].z())) driftFlux_[iCell].z() = -0.05*driftFlux_min[iCell].z();
	     }
	   else
	     {
	       if (driftFlux_[iCell].z() < (-0.05*driftFlux_min[iCell].z())) driftFlux_[iCell].z() = -0.05*driftFlux_min[iCell].z();
	       if (driftFlux_[iCell].z() > (0.93*driftFlux_min[iCell].z())) driftFlux_[iCell].z()=0.93*driftFlux_min[iCell].z();
	     }
       
       /*
     	   if (driftFlux_[iCell].z() < driftFlux_min[iCell].z()) driftFlux_[iCell].z() = driftFlux_min[iCell].z();
     	   if (driftFlux_[iCell].z() > 0) driftFlux_[iCell].z()=0;
       */
      }
      
  }
  
  driftFlux_.correctBoundaryConditions();

}

