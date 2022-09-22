/*---------------------------------------------------------------------------*\
particleSlipJohnsonJacksonFvPatchVectorField

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
    Foam::particleSlipJohnsonJacksonFvPatchVectorField

SourceFiles
    particleSlipJohnsonJacksonFvPatchVectorField.C
    
Description
    Johnson and Jackson boundary condition for the velocity of the particle
    phase.
    
Authors
    Alberto Passalacqua <albertop@iastate.edu; albert.passalacqua@gmail.com>
    Juho Peltola <juho.peltola@tut.fi>

Notes
    The code is known to work with OpenFOAM 1.7.x

\*---------------------------------------------------------------------------*/

#ifndef particleSlipJohnsonJacksonFvPatchVectorField_H
#define particleSlipJohnsonJacksonFvPatchVectorField_H

#include "fvPatchFields.H"
#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class particleSlipJohnsonJacksonFvPatch Declaration
\*---------------------------------------------------------------------------*/


class particleSlipJohnsonJacksonFvPatchVectorField
:
    public partialSlipFvPatchVectorField
{
    // Private data
	scalar specularityCoefficient_;	// Specularity coefficient

public:

    //- Runtime type information
    TypeName("particleSlipJohnsonJackson");


    // Constructors

        //- Construct from patch and internal field
        particleSlipJohnsonJacksonFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        particleSlipJohnsonJacksonFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  particleSlipJohnsonJacksonFvPatchVectorField onto a new patch
        particleSlipJohnsonJacksonFvPatchVectorField
        (
            const particleSlipJohnsonJacksonFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new particleSlipJohnsonJacksonFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        particleSlipJohnsonJacksonFvPatchVectorField
        (
            const particleSlipJohnsonJacksonFvPatchVectorField&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const DimensionedField<Vector<double>, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new particleSlipJohnsonJacksonFvPatchVectorField(*this, iF)
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
