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

#include "averageSlipVel.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

#define SMALL 1e-30

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(averageSlipVel, 0);

addToRunTimeSelectionTable
(
    forceModel,
    averageSlipVel,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
averageSlipVel::averageSlipVel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    fluidVelFieldName_(propsDict_.lookup("fluidVelFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (fluidVelFieldName_)),
    particleVelFieldName_(propsDict_.lookup("particleVelFieldName")),
    Us_(sm.mesh().lookupObject<volVectorField> (particleVelFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    rhoFluidName_(propsDict_.lookup("rhoFluidName")),
    rhoFluid_(sm.mesh().lookupObject<volScalarField> (rhoFluidName_)),
    rhoParticle_(readScalar(propsDict_.lookup("rhoParticle"))),
    outputDirName_(propsDict_.lookup("outputDirName"))
{
    //Initialize the output stream
    fileName probeDir = outputDirName_;
    probeDir = probeDir/U_.mesh().time().timeName();
    mkDir(probeDir);
  
    sPtrUSlip_ =  new OFstream(probeDir/"uSlip");
    *sPtrUSlip_  << "#Time    aveVoidfraction  uSlip.x uSlip.y uSlip.z" << endl;
    
    sPtrIntMomentum_ = new OFstream(probeDir/"integralMomentum");
    *sPtrIntMomentum_ << "#Time  refMomentum  totalMomentum.x totalMomentum.y totalMomentum.z" << endl;
    
    sPtrVelStats_    = new OFstream(probeDir/"velStats");
    *sPtrVelStats_    << "#Time  AveUFluid | AveUParticle" << endl;
    
    //Misc Warnings
    Warning << "implicit and explicit momentum Coupling must be used" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

averageSlipVel::~averageSlipVel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void averageSlipVel::setForce() const
{
    scalar ds(0);
    scalar Volp(0);
    vector Up(0,0,0);
    vector sumUp(0,0,0); 
    vector particleMomentum(0,0,0);
    scalar nP(0);

    for(int index = 0;index < particleCloud_.numberOfParticles(); index++)
    {
        ds   = 2*particleCloud_.radius(index);
        Volp = ds*ds*ds*M_PI/6;
        Up = particleCloud_.velocity(index);
        // Particle velocity
        sumUp += Up;
        // Particle momentum
        particleMomentum += rhoParticle_*Volp*Up;
    } 
    // Average particle velocity
    //vector aveParticleVel = sumUp/particleCloud_.numberOfParticles();
    nP = particleCloud_.numberOfParticles();
    // Parallel computations
    reduce(sumUp, sumOp<vector>());      
    reduce(nP, sumOp<scalar>()); 
    vector aveParticleVel = sumUp/nP;    
    
    // Total volume of domain
    scalar volDomain(0);
    forAll(voidfraction_,cellI)
    {
      volDomain += voidfraction_.mesh().V()[cellI];
    }
    // Parallel computations
    reduce(volDomain, sumOp<scalar>());  

    // Domain integrated fluid volume fraction
    scalar aveVoidfraction = fvc::domainIntegrate(voidfraction_).value();

    // Domain averaged fluid velocity
    vector domainAveFluidVel = fvc::domainIntegrate(voidfraction_*U_).value()/aveVoidfraction; 

    // Averaged slip velocity
    vector uSlip = domainAveFluidVel - aveParticleVel;

    // Domain averaged fluid momentum
    vector fluidMomentum = fvc::domainIntegrate(rhoFluid_*voidfraction_*U_).value();
    vector integralMomentum = fluidMomentum + particleMomentum; 
    scalar refMomentum  = aveVoidfraction*rhoParticle_*mag(uSlip); //reference momentum
    
    //Finalize and write to file
    if(Pstream::master()) //Write only if master
    {
    *sPtrUSlip_ << setw(IOstream::defaultPrecision() + 6) 
                << U_.mesh().time().value()     << "   " 
                << 1-aveVoidfraction/volDomain  << "   "
                << uSlip[0]                     << "   " 
                << uSlip[1]                     << "   "
                << uSlip[2]                     << "   "
                << endl;

    *sPtrIntMomentum_ << setw(IOstream::defaultPrecision() + 6) 
                << U_.mesh().time().value() << "   " 
                << refMomentum << "   "
                << integralMomentum[0] << "   "
                << integralMomentum[1] << "   "
                << integralMomentum[2] << "   "
                << endl;
                
    *sPtrVelStats_   << setw(IOstream::defaultPrecision() + 6) 
                << U_.mesh().time().value() << "   " 
                << domainAveFluidVel[0]<< "   "
                << domainAveFluidVel[1]<< "   "
                << domainAveFluidVel[2]<< "   "
                << "|  "
                << aveParticleVel[0]<< "   "
                << aveParticleVel[1]<< "   "
                << aveParticleVel[2]<< "   "
                << endl;               
    }
 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

