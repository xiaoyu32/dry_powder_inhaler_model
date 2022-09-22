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
#include "cfdemCloudStress.H"
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
cfdemCloudStress::cfdemCloudStress
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    velCorrection_(NULL),
    collStress_(NULL),
    sigmaCollCell_
    (
	IOobject
	(
            "sigmaColl",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	 ),
	mesh,
	dimensionedSymmTensor( "zero", dimensionSet(1,-1,-2,0,0), symmTensor(0,0,0,0,0,0) )
     )    
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudStress::~cfdemCloudStress()
{
    dataExchangeM().destroy(velCorrection_,3);
    dataExchangeM().destroy(collStress_,6);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
void Foam::cfdemCloudStress::giveDEMdata()
{
    //cfdemCloud::giveDEMdata();
    dataExchangeM().giveData("velocityCorrection","vector-atom",velCorrection_);
}
*/

bool Foam::cfdemCloudStress::reAllocArrays() const
{    
    if(cfdemCloud::reAllocArrays())
    {
        dataExchangeM().allocateArray(velCorrection_,0.,3);
        dataExchangeM().allocateArray(collStress_,0.,6);	
    }    
    return true;
}



void Foam::cfdemCloudStress::calcSendVelocityCorrection
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
    //cfdemCloudStress::reAllocArrays();
    for(int index=0; index< numberOfParticles(); index++)
    {
        for(int i=0;i<3;i++) velCorrection_[index][i] = -velocityCorrection[0][i];
    }    	
    dataExchangeM().giveData("velocity_correction","vector-atom",velCorrection_);
    
}

void Foam::cfdemCloudStress::getDEMStressData()
{
    dataExchangeM().getData("collStress","vector-atom",collStress_);		
}

symmTensor Foam::cfdemCloudStress::collStress(int index)
{
    symmTensor coll;
    for(int i=0;i<6;i++) coll[i] = collStress_[index][i];
    return coll;
}

void Foam::cfdemCloudStress::writeCollStress
(
	fvMesh& mesh,
	cfdemCloud& particleCloud_	
) 
{
    label cellID(-1);
    symmTensor coll(0,0,0,0,0,0);
    
    // Reset
    forAll(mesh.C(),cellI)
    {
	for(int ii=0; ii<6; ii++)
       	 sigmaCollCell_[cellI][ii] = 0.;		    	    
    }    

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
       cellID = particleCloud_.cellIDs()[index][0];
       coll = collStress(index);
 	
       if( cellID>-1 )
       {
        for(int ii=0; ii<6; ii++)
       	 sigmaCollCell_[cellID][ii] += coll[ii];			   
       }
    }
    
    forAll(mesh.C(),cellI)
    {
	for(int ii=0; ii<6; ii++)
       	 sigmaCollCell_[cellI][ii] /= mesh.V()[cellI];		    	    
    }
    
    sigmaCollCell_.write();

    
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
