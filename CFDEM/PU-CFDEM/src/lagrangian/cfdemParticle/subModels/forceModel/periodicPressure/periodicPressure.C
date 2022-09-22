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

#include "error.H"

#include "periodicPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(periodicPressure, 0);

addToRunTimeSelectionTable
(
    forceModel,
    periodicPressure,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
periodicPressure::periodicPressure
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    velocityFieldName_(propsDict_.lookup("velocityFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velocityFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),     
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),      
    rhoParticle_(readScalar(propsDict_.lookup("rhoParticle"))), 
    gravityFieldName_(propsDict_.lookup("gravityFieldName")),
    pressureMultiplier_( 1.0 ),
    #if defined(version21) || defined(version16ext)
        g_(sm.mesh().lookupObject<uniformDimensionedVectorField> (gravityFieldName_))
    #elif defined(version15)
        g_(dimensionedVector(sm.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(gravityFieldName_)).value())
    #endif 
{    
    if (propsDict_.found("verbose")) verbose_=true;
    
    if ( propsDict_.found("pressureMultiplier") )
    {
        pressureMultiplier_= readScalar( propsDict_.lookup("pressureMultiplier") );     
	Info << "Source term multiplier = " << pressureMultiplier_  << endl; 
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

periodicPressure::~periodicPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void periodicPressure::setForce() const
{
    scalar ds(0);
    scalar Volp(0);
    vector Up(0,0,0);
    scalar volParts(0);
    scalar massParts(0);
    vector momParts(0,0,0);
    label cellID(-1);   
 
    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        ds   = 2*particleCloud_.radius(index);
        Volp = ds*ds*ds*M_PI/6;
	Up = particleCloud_.velocity(index);

        cellID = particleCloud_.cellIDs()[index][0];
        // Some stupid periodic boundary problem 
        //if( cellID>-1 )
        {	
	  // Total volume of particles
	  volParts += Volp;
	  // Total mass of particles
	  massParts+= rhoParticle_*Volp;
        }
    }    	

    // Calculated weighted of fluid & particles
    scalar volDomain(0);
    scalar aveRhoFluid(0);
    scalar volFluid(0);
    vector totalFluidWeight (0,0,0);
    
    forAll(voidfraction_,cellI)
    {
        volDomain          += voidfraction_.mesh().V()[cellI];
	aveRhoFluid	   += voidfraction_.mesh().V()[cellI]*rho_[cellI];
	volFluid	   += voidfraction_.mesh().V()[cellI]*voidfraction_[cellI];
	totalFluidWeight   += voidfraction_.mesh().V()[cellI]*voidfraction_[cellI]*rho_[cellI] *g_.value();	
    }   
    
    // Parallel computations
    reduce(volDomain         , sumOp<scalar>());
    reduce(aveRhoFluid	     , sumOp<scalar>());
    reduce(volFluid	     , sumOp<scalar>());
    reduce(totalFluidWeight  , sumOp<vector>());    
    
    // Domain-averaged solid volume fraction in domain, total_vol_particle/vol_domain
    scalar domainPartVolFraction = volParts/volDomain;
    // Parallel computations
    reduce(domainPartVolFraction ,sumOp<scalar>());
    // Domain-averaged fluid volume fraction in domain
    scalar domainFluidVolFraction = 1. - domainPartVolFraction; //volFluid/volDomain; 
    // Domain-averaged density of fluid
    aveRhoFluid /= volDomain;   
        
    // Total weight in the domain
    volVectorField source_
    (   
        IOobject
        (
            "source",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	particleCloud_.mesh(),
        dimensionedVector("source", dimensionSet(1,-2,-2,0,0), pressureMultiplier_ * ( domainPartVolFraction*rhoParticle_*g_.value() + totalFluidWeight/volDomain ) ) 
    );
        
    if(verbose_)
    {
    	Info << "Total volume of domain = " << volDomain << endl;
	Info << "Domain-averaged solid volume fraction = " << domainPartVolFraction  << endl; 
	Info << "Domain-averaged fluid volume fraction = " << domainFluidVolFraction << endl;  	
	Info << "Source term = " << -fvc::domainIntegrate(source_).value()/volDomain  << endl; 
    }	

    //Set the Source in the Explicit Coupling Force Field
    particleCloud_.momCoupleM(1).setSourceField(source_); //set source in the explicit coupling force		
		
    label cellI=0; ds = 0; Volp = 0 ;
    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle found on this processor
            {
                //Calc the particle volume
                ds   = 2*particleCloud_.radius(index);
                Volp = ds*ds*ds*M_PI/6;

                // set force on particle
                for(int j=0;j<3;j++) 
                {
                    // calc particle's static pressure gradient force
		    DEMForces()[index][j] -= Volp*source_[cellI][j];
                }
            	if(verbose_ && index <1)  Info << "index = " << index << " Static pressure gradient force ("
                                                            << DEMForces()[index][0] << " "
                                                            << DEMForces()[index][1] << " "
                                                            << DEMForces()[index][2] << ")" << endl;


            }
        //}
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
