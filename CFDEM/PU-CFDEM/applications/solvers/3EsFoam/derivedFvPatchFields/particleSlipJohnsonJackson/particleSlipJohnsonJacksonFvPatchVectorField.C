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

\*---------------------------------------------------------------------------*/

#include "particleSlipJohnsonJacksonFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

particleSlipJohnsonJacksonFvPatchVectorField::particleSlipJohnsonJacksonFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    specularityCoefficient_(p.size())
{}


particleSlipJohnsonJacksonFvPatchVectorField::particleSlipJohnsonJacksonFvPatchVectorField
(
    const particleSlipJohnsonJacksonFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(tdpvf, p, iF, mapper),
    specularityCoefficient_(tdpvf.specularityCoefficient_)
{}


particleSlipJohnsonJacksonFvPatchVectorField::particleSlipJohnsonJacksonFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    specularityCoefficient_(readScalar(dict.lookup("specularityCoefficient")))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        partialSlipFvPatchVectorField::evaluate();
    }
}


particleSlipJohnsonJacksonFvPatchVectorField::particleSlipJohnsonJacksonFvPatchVectorField
(
    const particleSlipJohnsonJacksonFvPatchVectorField& tdpvf,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(tdpvf, iF),
    specularityCoefficient_(tdpvf.specularityCoefficient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void particleSlipJohnsonJacksonFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if ((specularityCoefficient_ < 0) || (specularityCoefficient_ > 1))
    {
	FatalErrorIn
        (
            "particleSlipJohnsonJacksonFvPatchScalarField::"
            "updateCoeffs()"
        )   << "The value of the specularity coefficient has to be between 0 and 1."
            << abort(FatalError);
    }
    
    const dictionary& transportProperties = db().lookupObject<IOdictionary>
    (
        "transportProperties"
    );

    const dictionary& kineticTheoryProperties = db().lookupObject<IOdictionary>
    (
        "kineticTheoryProperties"
    );

    dictionary phaseaDictionary
    (
        transportProperties.subDict("phasea")
    );

    dimensionedScalar rhoa(phaseaDictionary.lookup("rho"));
    dimensionedScalar alphaMax(kineticTheoryProperties.lookup("alphaMax"));

    const fvPatchScalarField& alpha = 
	patch().lookupPatchField<volScalarField, scalar>("alpha");

    const fvPatchScalarField& g0 = 
	patch().lookupPatchField<volScalarField, scalar>("gs0");

    const fvPatchScalarField& mua = 
	patch().lookupPatchField<volScalarField, scalar>("mua");

    scalarField alphaPatch = alpha.patchInternalField() + 1.0e-6;
    scalarField ThetaPatch = max(alphaPatch, 1.0e-6);
    scalarField g0Patch = g0.patchInternalField();
    scalarField muaPatch = mua.patchInternalField();
    
    if (db().foundObject<volScalarField>("Theta"))
    {
	const fvPatchScalarField& Theta = 
	  patch().lookupPatchField<volScalarField, scalar>("Theta");
	
	ThetaPatch = Theta.patchInternalField();
    }
    
    // The partial slip BC in OpenFOAM is implemented as
    //
    // valueFraction*U + (1-valueFraction)*grad(U) = 0
    //
    // To find valueFraction, we re-write Johnson and Jackson BC as
    //
    // c*U + grad(U) = 0
    //
    // where
    //
    // c = valueFraction/(1 - valueFraction)
    //
    // As a consequence
    //
    // valueFraction = 1/(c + 1)
    
    scalarField c = (6.0*muaPatch*alphaMax.value())/
		    (M_PI*rhoa.value()*alphaPatch*g0Patch*
		     specularityCoefficient_*sqrt(3.0*ThetaPatch));

    this->valueFraction() = scalar(1)/(c + scalar(1));

    partialSlipFvPatchVectorField::updateCoeffs();
}

// Write
void particleSlipJohnsonJacksonFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("specularityCoefficient")
        << specularityCoefficient_ << token::END_STATEMENT << nl;
    writeEntry("value", os);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    particleSlipJohnsonJacksonFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
