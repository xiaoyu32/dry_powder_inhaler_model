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

#include "linearViscousStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameWithName(linearViscousStress, "linearViscousStress");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearViscousStress::linearViscousStress
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName,
    const volScalarField rho
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),

    ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce",
            coeffDict_,
            1.048
        )
    ),

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
//    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//++++++++++++++++++++++++++++++++++++++++++++++++++

tmp<fvVectorMatrix>linearViscousStress::divDevRhoReff
(
    const volScalarField& rho,
    const volScalarField alpha,
    volVectorField& U
) const
{
    return
    (
      - fvc::div((alpha*rho*nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(alpha*rho*nuEff(), U)
    );
}
//++++++++++++++++++++++++++++++++++++++++++++++++++
tmp<volSymmTensorField>linearViscousStress::devRhoReff(    const volScalarField& rho,
    const volScalarField alpha,
    volVectorField& U
) const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
       
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (
             "devRhoReff", (-alpha*rho*nuEff())
           *dev(twoSymm(fvc::grad(U)))
            )
	)
  );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++

void linearViscousStress::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct(gradU);
}


bool linearViscousStress::read()
{
    if (LESModel::read())
    {
        ce_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
