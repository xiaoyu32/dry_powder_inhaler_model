/*---------------------------------------------------------------------------*\
particleThetaJohnsonJacksonFvPatchScalarField

Copyright Information
    Copyright (C) 2008 Alberto Passalacqua 
    Copyright (C) 2008-2010 Juho Peltola
    
License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::particleThetaJohnsonJacksonFvPatchScalarField

SourceFiles
    particleThetaJohnsonJacksonFvPatchScalarField.C
    
Description
    Johnson and Jackson boundary condition for the granular temperature.
    
Authors
    Alberto Passalacqua <albertop@iastate.edu; albert.passalacqua@gmail.com>
    Juho Peltola <juho.peltola@tut.fi>

Notes
    The code is known to work with OpenFOAM 1.7.x

\*---------------------------------------------------------------------------*/

#ifndef particleThetaJohnsonJacksonFvPatchScalarField_H
#define particleThetaJohnsonJacksonFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class particleThetaJohnsonJacksonFvPatch Declaration
\*---------------------------------------------------------------------------*/

class particleThetaJohnsonJacksonFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data
    
        scalar restitutionCoefficient_;	// Particle-wall restitution coefficient
	scalar specularityCoefficient_;	// Specularity coefficient


public:

    //- Runtime type information
    TypeName("particleThetaJohnsonJackson");


    // Constructors

        //- Construct from patch and internal field
        particleThetaJohnsonJacksonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<double, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        particleThetaJohnsonJacksonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<double, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  particleThetaJohnsonJacksonFvPatchScalarField
        //  onto a new patch
        particleThetaJohnsonJacksonFvPatchScalarField
        (
            const particleThetaJohnsonJacksonFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<double, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new particleThetaJohnsonJacksonFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        particleThetaJohnsonJacksonFvPatchScalarField
        (
            const particleThetaJohnsonJacksonFvPatchScalarField&,
            const DimensionedField<double, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(const DimensionedField<double, volMesh>& iF) const
        {
            return tmp<fvPatchScalarField>
            (
                new particleThetaJohnsonJacksonFvPatchScalarField(*this, iF)
            );
        }



    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


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
