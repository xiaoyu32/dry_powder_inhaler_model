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

Class
    Foam::electrostaticsKTModel

SourceFiles
    electrostaticsKTModel.C

\*---------------------------------------------------------------------------*/

#ifndef electrostaticsKTModel_H
#define electrostaticsKTModel_H


#include "kineticTheoryModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class electrostaticsKTModel Declaration
\*---------------------------------------------------------------------------*/

class electrostaticsKTModel
{
    // Private data

        const phaseModel& phase1_;
        const volVectorField& U1_;
        const volScalarField& alpha1_;

        const dimensionedScalar& rho1_;
        const dimensionedScalar& da_;

        //- electostaticsKT dictionary 
        IOdictionary electrostaticsKTProperties_;

        //- use electostaticsKT or not.
        Switch electrostaticsKT_;

        //- Vacumm permittivity
        const dimensionedScalar epsilon0_;
	
	//- Young modulus
	const dimensionedScalar YoungModulus_;

	//- Poisson ratio
	const dimensionedScalar PoissonRatio_;
		
        //- KT dictionary 
        IOdictionary kineticTheoryProperties_;

        //- maximum packing
        const dimensionedScalar alphaMax_;

        //- Radial function
        autoPtr<kineticTheoryModels::radialModel> radialModel_;

        //- The radial distribution function
        volScalarField gs0E_;
	
	//- The granular energy/temperature
        volScalarField Theta_;

        //- The tribocharging conductivity
        volScalarField sigmaq_;

        //- The tribocharging diffusivity
        volScalarField kappaq_;

        //- Turbulent Schmidt Number  
        Switch chargeTurbulentDiffusivity_;
	
        //- The tribocharging turbulent diffusivity
        volScalarField kappaqt_;
	
	//- Phase fluctuation viscosity
	dimensionedScalar mupt_;
	
        //- Phase diffusion based on Theta,alpha,dp
        volScalarField DTheta_;	

        //- Radial function
        autoPtr<kineticTheoryModels::viscosityModel> viscosityModel_;
	
        //- coefficient of restitution
        const dimensionedScalar e_;
	
        //- Phase viscosity
        volScalarField mu1E_;
	
	//- theta (Jari's fav)
	volVectorField sigmaqEf_;	

    // Private Member Functions

        //- Disallow default bitwise copy construct
        electrostaticsKTModel(const electrostaticsKTModel&);

        //- Disallow default bitwise assignment
        void operator=(const electrostaticsKTModel&);


public:

    // Constructors

        //- Construct from components
        electrostaticsKTModel
        (
            const phaseModel& phase1,
            const volScalarField& alpha1
        );

    //- Destructor
    virtual ~electrostaticsKTModel();


    // Member Functions

        void compute();

        bool on() const
        {
            return electrostaticsKT_;
        }

        const volScalarField& sigmaq() const
        {
            return sigmaq_;
        }

        const volScalarField& kappaq() const
        {
            return kappaq_;
        }

        bool onTurbulentDiffusivity() const
        {
            return chargeTurbulentDiffusivity_;
        }

        const volScalarField& kappaqt() const
        {
            return kappaqt_;
        }

        const volVectorField& sigmaqEf() const
        {
            return sigmaqEf_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
