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
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(forceModel, 0);

defineRunTimeSelectionTable(forceModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forceModel::forceModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    treatExplicit_(false),
    treatDEM_(false),
    implDEM_(false),
    impParticleForces_
    (   IOobject
        (
            "impParticleForces",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,1,-2,0,0), vector(0,0,0)) // N
    ),
    expParticleForces_
    (   IOobject
        (
            "expParticleForces",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,1,-2,0,0), vector(0,0,0)) // N
    ),
    coupleForce_(true),
    modelType_(sm.modelType()),
//_AO
    probeIt_(sm.probeM().active()) //,
//_AO 06/02/2014 
/* 
//    cg_(1.)
*/
//_AO 
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceModel::~forceModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
/*tmp<volScalarField> forceModel::provideScalarField()
{
Info << "now providing a scalar field" << endl;
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "xxx",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimensionSet(0,0,0,0,0),
                0.0
            )
        )
    );

    manipulateScalarField(tsource());
    return tsource;
};*/

void forceModel::manipulateScalarField(volScalarField& field) const
{
    Info << "no scalar manipulation done" << endl;
    // do nothing
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
