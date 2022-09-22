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
#include "cfdemCloudPerio.H"
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

// Construct from components
cfdemCloudPerio::cfdemCloudPerio
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    velCorrection_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudPerio::~cfdemCloudPerio()
{
    dataExchangeM().destroy(velCorrection_,3);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
void Foam::cfdemCloudPerio::giveDEMdata()
{
    //cfdemCloud::giveDEMdata();
    dataExchangeM().giveData("velocityCorrection","vector-atom",velCorrection_);
}
*/

bool Foam::cfdemCloudPerio::reAllocArrays() const
{    
    if(cfdemCloud::reAllocArrays())
    {
        dataExchangeM().allocateArray(velCorrection_,0.,3);
    }    
    return true;
}

void Foam::cfdemCloudPerio::calcSendVelocityCorrection
(
    scalar& rhoParticle,
    volScalarField& rho,
    volScalarField& voidfraction,
    volVectorField& U,
    volVectorField& Us
)
{
    vector fluidMomentum = fvc::domainIntegrate(rho*voidfraction*U).value();
    vector solidMomentum = fvc::domainIntegrate(rhoParticle*(1-voidfraction)*Us).value();

    volVectorField velocityCorrection
    (   
	    IOobject
	    (
		"velCorr",
		U.mesh().time().timeName(),
		U.mesh(),
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedVector( "zero", dimensionSet(0,1,-1,0,0), 
	      ( fluidMomentum + solidMomentum ) 
	    / ( fvc::domainIntegrate(rhoParticle*(1-voidfraction)).value()+fvc::domainIntegrate(rho*voidfraction).value() ) )
    );

    Info << "Mixture velocity correction " << endl;
    Info << "Velocity correction at cell[0] = " << velocityCorrection[0] << endl;
    Info << " " << endl;

    // Total momentum correction for fluid and solid
    U = U-velocityCorrection;
    
    // Send to DEM part
    //cfdemCloudPerio::reAllocArrays();
    for(int index=0; index< numberOfParticles(); index++)
    {
        for(int i=0;i<3;i++) velCorrection_[index][i] = -velocityCorrection[0][i];
    }    	
    dataExchangeM().giveData("velocity_correction","vector-atom",velCorrection_);
    
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
