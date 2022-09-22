/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "noFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noFrictionalStress, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        noFrictionalStress,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noFrictionalStress::noFrictionalStress
(
    const dictionary& dict
)
:
    frictionalStressModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::noFrictionalStress::~noFrictionalStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::noFrictionalStress::
frictionalPressure
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{

    return
        scalar(0.0)*
        Fr*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta)
       /pow(max(alphaMax - alpha1, scalar(5.0e-2)), p);
}


Foam::tmp<Foam::volScalarField> Foam::noFrictionalStress::
frictionalPressurePrime
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{
    return 
	scalar(0.0)*
	Fr*
    (
        eta*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta - 1.0)
       *(alphaMax-alpha1)
      + p*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta)
    )/pow(max(alphaMax - alpha1, scalar(5.0e-2)), p + 1.0);
}


Foam::tmp<Foam::volScalarField> Foam::noFrictionalStress::muf
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D,
    const dimensionedScalar& phi
) const
{
    return scalar(0.0)* 
	dimensionedScalar("0.5", dimTime, 0.5)*pf*sin(phi);
}


// ************************************************************************* //
