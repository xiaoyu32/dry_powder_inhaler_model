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
#include "cfdemCloudScalar.H"
#include "voidFractionModel.H"
#include "forceModel.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
CfdemCloudScalar::CfdemCloudScalar
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh)
{ }
*/

// Construct from components
CfdemCloudScalar::CfdemCloudScalar
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    passiveScalarDensityModel_
    (
	passiveScalarDensityModel::New
	(
	   couplingProperties_,
	   *this
	)
    ),
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
     ),
     pscalar_( passiveScalarDensityM().numberOfPassiveScalars() ),
     dpscalar_( passiveScalarDensityM().numberOfPassiveScalars() ),
     numberOfPassiveScalars_( passiveScalarDensityM().numberOfPassiveScalars() )   
{ 

    for( int ii = 0; ii < numberOfPassiveScalars_; ++ii )
    {
        pscalar_[ii] = NULL;
	dpscalar_[ii] = NULL;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

CfdemCloudScalar::~CfdemCloudScalar()
{

  
    for( int ii = 0; ii < numberOfPassiveScalars_; ++ii )
    {
       dataExchangeM().destroy(pscalar_[ii], 1);
       dataExchangeM().destroy(dpscalar_[ii], 1);  
    }
    
    dataExchangeM().destroy(velCorrection_,3);
    dataExchangeM().destroy(collStress_,6);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::CfdemCloudScalar::getDEMdata()
{
    cfdemCloud::getDEMdata();
}

bool Foam::CfdemCloudScalar::reAllocArrays() const
{
    
    if(cfdemCloud::reAllocArrays())
    {
    
        //numberOfPassiveScalars_ = 2;
    
        Info<<"Reallocating "<<numberOfPassiveScalars_<<" scalars.."<<endl;

	for( int ii = 0; ii < numberOfPassiveScalars_; ++ii )
	{
	    dataExchangeM().allocateArray( pscalar_[ii],0., 1 );	    
            dataExchangeM().allocateArray( dpscalar_[ii],0., 1 );
	}	
	
	//dataExchangeM().allocateArray(potential_ ,0.,voidFractionM().maxCellsPerParticle());
	// For periodic box
        dataExchangeM().allocateArray(velCorrection_,0.,3);
        dataExchangeM().allocateArray(collStress_,0.,6);	

    }
    
    return true;
}

passiveScalarDensityModel& Foam::CfdemCloudScalar::passiveScalarDensityM()
{
    return passiveScalarDensityModel_();
}



//bool Foam::CfdemCloudScalar::evolve(volScalarField& rhoe)
void Foam::CfdemCloudScalar::evolvePassiveScalar
(
	List< autoPtr<volScalarField> >& passiveScalar,
	List< autoPtr<volScalarField> >& expPassiveScalarSource,
	List< autoPtr<volScalarField> >& impPassiveScalarSource,
	volScalarField& alpha,
	volVectorField& Us,
	volVectorField& U
)
{
    //numberOfParticlesChanged_ = false;
    //arraysReallocated_=false;
    
    //reAllocArrays();
    
    
    if ( cfdemCloud::evolve(alpha,Us,U) )
    {    	
    
	Info<<"Coupling passive scalar ( "<<numberOfPassiveScalars_<<" scalars )...\n"<<endl;
	
        
	for( int ii = 0; ii < numberOfPassiveScalars_; ++ii )
	{
	    word iiNum( name( ii ) );	    
	    word scalarName( "passiveScalar" + iiNum );
	     
            dataExchangeM().getData( scalarName,"scalar-atom",pscalar_[ii]);
	}
		
	passiveScalarDensityM().getScalarTransfer( 
							pscalar_, 
							expPassiveScalarSource, 
							impPassiveScalarSource,
							alpha,
							Us,
							U 
						 );
		
    }
        
    Info<<"Coupling Done."<<endl;
    
    // if(verbose_)    #include "debugInfo.H";

}

void Foam::CfdemCloudScalar::updatePassiveScalarSource
( 
	List< autoPtr<volScalarField> >& passiveScalar
)
{
    
    //send passive scalar flux to particles
    Info<<"Sending passive scalar flux to particles ( "<<numberOfPassiveScalars_<<" scalars )...\n"<<endl;
    
    passiveScalarDensityM().computeScalarTransfer( pscalar_, dpscalar_, passiveScalar );
    
    for( int ii = 0; ii < numberOfPassiveScalars_; ++ii )
    {
	word iiNum( name( ii ) );	    
	word scalarFluxName( "passiveScalarFlux" + iiNum );
	
	dataExchangeM().giveData( scalarFluxName,"scalar-atom", dpscalar_[ii] );    
    }
    	
    // do particle IO
    IOM().dumpDEMdata();

}


//- Calculate correction and send it to LIGGGHTS

void Foam::CfdemCloudScalar::calcSendVelocityCorrection
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
    //CfdemCloudScalar::reAllocArrays();
    for(int index=0; index< numberOfParticles(); index++)
    {
        for(int i=0;i<3;i++) velCorrection_[index][i] = -velocityCorrection[0][i];
    }    	
    dataExchangeM().giveData("velocity_correction","vector-atom",velCorrection_);
    
}

void Foam::CfdemCloudScalar::getDEMStressData()
{
    dataExchangeM().getData("collStress","vector-atom",collStress_);		
}

symmTensor Foam::CfdemCloudScalar::collStress(int index)
{
    symmTensor coll;
    for(int i=0;i<6;i++) coll[i] = collStress_[index][i];
    return coll;
}

void Foam::CfdemCloudScalar::writeCollStress
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
