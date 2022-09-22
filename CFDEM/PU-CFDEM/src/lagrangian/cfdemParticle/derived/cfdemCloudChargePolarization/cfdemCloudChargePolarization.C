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
#include "chargeDensityModel.H"
#include "cfdemCloudChargePolarization.H"
#include "voidFractionModel.H"
#include "forceModel.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "polarizationModel.H"
#include "permittivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
CfdemCloudChargePolarization::CfdemCloudChargePolarization
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh)
{ }
*/

// Construct from components
CfdemCloudChargePolarization::CfdemCloudChargePolarization
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    chargeDensityModel
    (
	chargeDensityModel::New
	(
	   couplingProperties_,
	   *this
	)
    ),
    polarizationModel
    (
	PolarizationModel::New
	(
	   couplingProperties_,
	   *this
	)
    ),
    permittivityModel
    (
	PermittivityModel::New
	(
	   couplingProperties_,
	   *this
	)    
    ),
    surfaceChargeDensity_
    (
         IOobject
        (
            "surfaceChargeDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedScalar( "zero", dimensionSet(0,-2,0,0,0),scalar(0) )      
    ),
    wallDataExchangeModel_(),
    wallSurfaceChargeDensity_( "surfaceChargeDensity" ), 
    normalEfieldContainer_( "normalElectricField" ),   
    pcharge_(NULL),
    ef_grad(NULL),
    pef_(NULL),
    insulatingWallsFlag_
    (
       couplingProperties_.lookupOrDefault<Switch>("insulator", false )
    ) 
{ 
    
    // -- if simulation has insulating surfaces, initialize the wall data transfer model --
    if( insulatingWallsFlag_ )
    {
        wallDataExchangeModel_.set
	( 
	   wallDataExchangeModel::New
	   (
	      dataExchangeM().getLmp(),
	      *this
	   ).ptr()	
	);
    }
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

CfdemCloudChargePolarization::~CfdemCloudChargePolarization()
{
    dataExchangeM().destroy(pcharge_,1);
    dataExchangeM().destroy(pef_,3);
    dataExchangeM().destroy(ef_grad,9);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::CfdemCloudChargePolarization::getDEMdata()
{
    cfdemCloud::getDEMdata();
}

bool Foam::CfdemCloudChargePolarization::reAllocArrays() const
{
    
    if(cfdemCloud::reAllocArrays())
    {
        dataExchangeM().allocateArray(pcharge_,0.,1);
        dataExchangeM().allocateArray(pef_,0.,3);
	// For periodic box	
        dataExchangeM().allocateArray(ef_grad,0.,9);
	
	return true;
    }
    
    return false;
}

chargeDensityModel& Foam::CfdemCloudChargePolarization::chargeDensityM()
{
    return chargeDensityModel();
}

PolarizationModel& Foam::CfdemCloudChargePolarization::polarizationM()
{
    return polarizationModel();
}

PermittivityModel& Foam::CfdemCloudChargePolarization::permittivityM()
{
    return permittivityModel();
}


//bool Foam::CfdemCloudChargePolarization::evolve(volScalarField& rhoe)
void Foam::CfdemCloudChargePolarization::evolve
(
   volVectorField& Ef,
   volScalarField& rhoe,
   volScalarField& alpha,
   volVectorField& Us,
   volVectorField& U,
   volScalarField& permittivity
)
{
    //numberOfParticlesChanged_ = false;
    //arraysReallocated_=false;
    
    if ( cfdemCloud::evolve(alpha,Us,U) )
    {    	
	Info<<"Coupling Charge...\n"<<endl;
	
        dataExchangeM().getData("charge","scalar-atom",pcharge_);

	chargeDensityM().resetChargeDensity();
    	chargeDensityM().setChargeDensity( pcharge_ );
	permittivityM().setPermittivity( alpha, permittivity ); // -- update permittivity --	
	
	// -- get boundary charge density --
	if( insulatingWallsFlag_ )
	{
	    // -- get DEM wall surface charge density --
	    wallDataExchangeM().getData( wallSurfaceChargeDensity_ );

	    // -- assign DEM surface charge density to the corresponding CFD surface field --
	    wallDataExchangeM().transferDataToField( wallSurfaceChargeDensity_, surfaceChargeDensity_ );	
	        
	}
	
    }
    
    
    rhoe.internalField() = chargeDensityM().chargeDensityInterp();	
    
    Info<<"Coupling Done."<<endl;
    
    // if(verbose_)    #include "debugInfo.H";

    // do particle IO
    IOM().dumpDEMdata();

}

void Foam::CfdemCloudChargePolarization::evolveElectricField( volVectorField& Ef, volTensorField& gradEf )
{

    Info<<"Coupling Electric Field..."<<endl;
    
    chargeDensityM().getElectricField( pef_, Ef );
    polarizationM().getElectricFieldGradient( ef_grad, gradEf );
    
    dataExchangeM().giveData("electricfield","vector-atom", pef_);
    dataExchangeM().giveData("electricfield_gradient","vector-atom", ef_grad);
    
    // -- compute normal electric field at surfaces & communicate to DEM walls for charge transfer --
    if( insulatingWallsFlag_ )
    {
       tmp<surfaceScalarField> normalEfield = ( fvc::interpolate( Ef ) & Ef.mesh().Sf() )/ Ef.mesh().magSf();
       
       // -- allocate memory on container for data transfer --
       wallDataExchangeM().initData( normalEfieldContainer_ );
       
       wallDataExchangeM().transferDataFromField( normalEfieldContainer_, normalEfield() );
       
       // -- send electric field information to DEM wall elements --
       wallDataExchangeM().giveData( normalEfieldContainer_ );
    }
    
    Info<<"Electricfield Coupled."<<endl;
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
