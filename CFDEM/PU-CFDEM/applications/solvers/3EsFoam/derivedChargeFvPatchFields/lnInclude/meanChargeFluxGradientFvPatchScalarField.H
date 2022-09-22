/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::meanChargeFluxFvPatchScalarField

Description
    gradient boundary condition for mean charge

SourceFiles
    meanChargeFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef meanChargeFluxFvPatchScalarField_H
#define meanChargeFluxFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class meanChargeFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class meanChargeFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the electric field
        word EfName_;

        //- Name of the granulat temperature field
        word ThetaName_;
	
	//- Name of phase-fraction field
        word phaseName_;
	
	//- Name of charge field
        word chargeName_;
	
	//- Work function difference with wall and solid phase
	scalar deltaWorkFunction_;
	
	//- Wall normal restitution coefficient
	scalar ew_;
	

public:

    //- Runtime type information
    TypeName("meanChargeFlux");


    // Constructors

        //- Construct from patch and internal field
        meanChargeFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        meanChargeFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given meanChargeFluxFvPatchScalarField
        // onto a new patch
        meanChargeFluxFvPatchScalarField
        (
            const meanChargeFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        meanChargeFluxFvPatchScalarField
        (
            const meanChargeFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new meanChargeFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        meanChargeFluxFvPatchScalarField
        (
            const meanChargeFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new meanChargeFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
