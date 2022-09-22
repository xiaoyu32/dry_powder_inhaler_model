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

#include "fileName.H"
#include "cfdemCloud.H"

#include "apiTransferModel.H"
#include "apiDensityModel.H"

#include "cfdemCoarseDPICloud.H"

#include "wallDataExchangeModel.H"
#include "voidFractionModel.H"
#include "forceModel.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
CfdemCloudCharge::CfdemCloudCharge
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh)
{ }
*/

// Construct from components
cfdemCoarseDPICloud::cfdemCoarseDPICloud
(
    const fvMesh& mesh,
    const volScalarField& alpha,
    const volVectorField& Us,
    const volVectorField& U
)
:
    cfdemCloud(mesh),
    carrierConcentration_(NULL),
    carrierFlux_(NULL),
    apiFluidDensity_
    (
        IOobject
        (
            "apiDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    apiCarrierDensity_
    (
        IOobject
        (
            "apiCarrierDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedScalar( "apiCarrierDensity_", dimensionSet(0,-3,0,0,0),scalar(0) )     
    ),
    apiImplSource_
    (
        IOobject
        (
            "apiImplSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh,
        dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0),scalar(0) )      
    ),
    apiExplSource_
    (
        IOobject
        (
            "apiExplSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh,
        dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0),scalar(0) )          
    ),
    apiWallDensity_
    (
         IOobject
        (
            "apiWallDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedScalar( "zero", dimensionSet(0,-2,0,0,0),scalar(0) )   
    ),
    apiDiffusivity_
    (
        IOobject
        (
            "apiDiffusivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedScalar( "one", dimensionSet(0,2,-1,0,0),scalar(1.0) )      
    ),
    alpha_( alpha ),
    Us_( Us ),
    U_( U ),
    demWallConcentration_( "apiSurfaceDensity" ),
    apiTransferModel_
    (
	apiTransferModel::New
	(
	   this->couplingProperties_.subDict( "apiTransferModelProps" ),
	   *this
	)    
    ),
    apiDensityModel_
    (
        apiDensityModel::New
	(
	   this->couplingProperties_,
	   *this
	) 
    ),
    wallDataExchangeModel_
    (
	wallDataExchangeModel::New
	(
	   dataExchangeM().getLmp(),
	   *this
	)
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCoarseDPICloud::~cfdemCoarseDPICloud()
{
    dataExchangeM().destroy(carrierConcentration_,1);
    dataExchangeM().destroy(carrierFlux_,1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cfdemCoarseDPICloud::getDEMdata()
{
    cfdemCloud::getDEMdata();  
}


bool Foam::cfdemCoarseDPICloud::reAllocArrays() const
{
    
    if(cfdemCloud::reAllocArrays())
    {
        dataExchangeM().allocateArray(carrierConcentration_,0.,1);
	dataExchangeM().allocateArray(carrierFlux_,0.,1);
	return true;
    }
    
    return false;
}


//bool Foam::CfdemCloudCharge::evolve(volScalarField& rhoe)
bool Foam::cfdemCoarseDPICloud::evolve
(
   volScalarField& alpha,
   volVectorField& Us,
   volVectorField& U
)
{
    
    bool returnValue = false;
    
    if ( cfdemCloud::evolve(alpha,Us,U) )
    {    
    
	Info<<"Coupling API particle concentration...\n"<<endl;
	// -- get API particle concentration and collisional flux --
	dataExchangeM().getData("carrierConcentration","scalar-atom",carrierConcentration_);
	dataExchangeM().getData("carrierCollisionFlux","scalar-atom",carrierFlux_);  

	apiCarrierDensity_.internalField() = 0;
	apiDensityM().apiDensityInterp( carrierConcentration_, apiCarrierDensity_ );
	
	// -- flux is given in terms of particles per unit volume ---> divide by the time step --
	//
	Info<< "Hello from ="<<" cfdemCoarseDPICloud.C"<<endl;
	apiExplSource_.internalField() = 0;
	apiDensityM().apiDensityInterp( carrierFlux_, apiExplSource_ );
	apiExplSource_.internalField() /= this->mesh().time().deltaT().value();
	
	// -- get DEM wall API particle wall surface density --
	wallDataExchangeM().getData( demWallConcentration_ );
	
	// -- assign DEM concentration to the CFD surface field --
	wallDataExchangeM().transferDataToField( demWallConcentration_, apiWallDensity_ );
	
	returnValue = true;
	
    }
        
    Info<<"Coupling Done."<<endl;
    
    return returnValue;
    
}

void Foam::cfdemCoarseDPICloud::updateExplSources()
{
    apiTransferM().explicitSources( apiExplSource_ );
}

void Foam::cfdemCoarseDPICloud::updateImplSources()
{
    apiImplSource_.internalField() = 0;
    apiTransferM().implicitSources( apiImplSource_ );
}

void Foam::cfdemCoarseDPICloud::commValues()
{
    
    cfdemCloud::commValues();
    
    // -- evolve API density field --
    this->evolveAPIfields();
    
    
    // -- interpolate wall surface concentration back to the DEM walls --
    wallDataExchangeM().transferDataFromField( demWallConcentration_, apiWallDensity_ );
    
    // -- send wall surface density to DEM --
    wallDataExchangeM().giveData( demWallConcentration_ );
    
    
    // -- send API particle flux to carrier particles to DEM --
    dataExchangeM().giveData("carrierFluidFlux","scalar-atom",carrierFlux_);
    
}

// -- evolve fluid API particle density fields --
void Foam::cfdemCoarseDPICloud::evolveAPIfields()
{
    
    Info<<"Solving API particle concentration..."<<endl;
    
    surfaceScalarField uflux = fvc::interpolate( alpha_ ) * linearInterpolate( U_ ) & mesh_.Sf();
    
    
    // -- update implicit and explicit sources --
    updateExplSources();
    updateImplSources();
    
    
    // -- update diffusivity (use turbulence diffusivity) --
    apiDiffusivity_ = this->turbulence().nut();
    
    fvScalarMatrix apiScalarEqn    
    (
        fvm::ddt( alpha_, apiFluidDensity_ ) 
      + fvm::div( uflux, apiFluidDensity_ )
        ==
      	fvm::laplacian( alpha_ * apiDiffusivity_, apiFluidDensity_ )
      - fvm::Sp( apiImplSource_, apiFluidDensity_ ) // -- implicit flux --
      + apiExplSource_ // -- explicit flux --
    );
    
    apiScalarEqn.relax();
    apiScalarEqn.solve();
    
    // -- compute API carrier flux from the fluid to carrier particles --
    this->resetArray(carrierFlux_,numberOfParticles(),1);
    apiTransferM().computeSource( apiFluidDensity_, carrierFlux_ );
    
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
