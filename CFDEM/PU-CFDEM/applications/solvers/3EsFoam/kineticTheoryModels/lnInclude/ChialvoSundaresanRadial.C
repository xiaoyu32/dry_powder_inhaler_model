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

#include "ChialvoSundaresanRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(ChialvoSundaresan, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        ChialvoSundaresan,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::ChialvoSundaresan
(
    const dictionary& dict
)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::~ChialvoSundaresan()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::g0
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{

    const scalar alpha2 = 0.58;
    const scalar alpha_c = alphaMax.value();
    const scalar alphastar = alpha_c - 0.01;
    return  
    	(1.0 - min(alphastar,alpha)/2.) / (pow(1.0 - min(alphastar,alpha), 3.)) 
		     + alpha2 * min(alphastar,alpha) * min(alphastar,alpha) / pow( alpha_c - min(alphastar,alpha) , 1.5) 		    
		     + (- 1.0 / ( 2.0 * pow(1.0 - alphastar, 3.))
             + 3.0 * ( 1.0 - alphastar/2.) / (pow(1.0 - alphastar, 4)) + 2.0 * alpha2 * alphastar / pow( alpha_c - alphastar ,1.5)
             + 3.0 / 2.0 * alpha2 * alphastar * alphastar / pow( alpha_c - alphastar, 2.5))*max(alpha-alphastar,0.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::g0prime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{
    const scalar alpha2 = 0.58;
    const scalar alpha_c = alphaMax.value();
    const scalar alphastar = alpha_c - 0.01;
    
    volScalarField valg0CSprime = 
    - 1.0 / ( 2.0 * pow(1.0 - min(alphastar,alpha), 3))
    + 3.0 * ( 1.0 - min(alphastar,alpha)/2.) / (pow(1.0 - min(alphastar,alpha), 4));
    
    return  
    	valg0CSprime  
    + 2.0 * alpha2 * min(alphastar,alpha) / pow( alpha_c - min(alphastar,alpha) ,1.5)
    + 3.0 / 2.0 * alpha2 * min(alphastar,alpha) * min(alphastar,alpha) / pow( alpha_c - min(alphastar,alpha), 2.5);

}


// ************************************************************************* //
