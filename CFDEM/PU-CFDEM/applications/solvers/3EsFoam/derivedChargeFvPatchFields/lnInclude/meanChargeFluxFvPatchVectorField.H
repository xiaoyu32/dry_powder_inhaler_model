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
    Foam::meanChargeFluxFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Charge flux boundary condition

\*---------------------------------------------------------------------------*/

#ifndef meanChargeFluxFvPatchVectorField_H
#define meanChargeFluxFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class meanChargeFluxFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class meanChargeFluxFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of the electric field
        word EfName_;

        //- Effect of electric field on charge transfer at the wall
        Switch EfBool_;

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
        meanChargeFluxFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        meanChargeFluxFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given meanChargeFluxFvPatchVectorField
        //  onto a new patch
        meanChargeFluxFvPatchVectorField
        (
            const meanChargeFluxFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        meanChargeFluxFvPatchVectorField
        (
            const meanChargeFluxFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new meanChargeFluxFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        meanChargeFluxFvPatchVectorField
        (
            const meanChargeFluxFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new meanChargeFluxFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
