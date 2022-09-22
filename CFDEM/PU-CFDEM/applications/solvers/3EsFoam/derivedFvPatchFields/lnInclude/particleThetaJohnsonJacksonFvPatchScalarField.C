/*---------------------------------------------------------------------------*\
particleThetaJohnsonJacksonFvPatchScalarField

Copyright Information
    Copyright (C) 2008 Alberto Passalacqua 
    Copyright (C) 2008-2010 Juho Peltola
    Copyright (C) 2010 Alberto Passalacqua 
    
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

#include "particleThetaJohnsonJacksonFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

particleThetaJohnsonJacksonFvPatchScalarField::particleThetaJohnsonJacksonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_(p.size()),
    specularityCoefficient_(p.size())
{
    this->refValue() = *this;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


particleThetaJohnsonJacksonFvPatchScalarField::particleThetaJohnsonJacksonFvPatchScalarField
(
    const particleThetaJohnsonJacksonFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


particleThetaJohnsonJacksonFvPatchScalarField::particleThetaJohnsonJacksonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_(readScalar(dict.lookup("restitutionCoefficient"))),
    specularityCoefficient_(readScalar(dict.lookup("specularityCoefficient")))
{
    this->refValue() = *this;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        particleThetaJohnsonJacksonFvPatchScalarField::evaluate();
    }
}


particleThetaJohnsonJacksonFvPatchScalarField::particleThetaJohnsonJacksonFvPatchScalarField
(
    const particleThetaJohnsonJacksonFvPatchScalarField& pivpvf,
    const DimensionedField<double, volMesh>& iF
)
:
    mixedFvPatchScalarField(pivpvf, iF),
    restitutionCoefficient_(pivpvf.restitutionCoefficient_),
    specularityCoefficient_(pivpvf.specularityCoefficient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleThetaJohnsonJacksonFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void particleThetaJohnsonJacksonFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void particleThetaJohnsonJacksonFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if ((restitutionCoefficient_ < 0) || (restitutionCoefficient_ > 1))
    {
	FatalErrorIn
        (
            "particleThetaJohnsonJacksonFvPatchScalarField::"
            "updateCoeffs()"
        )   << "The value of the restitution coefficient has to be between 0 and 1."
            << abort(FatalError);
    }

    if ((specularityCoefficient_ < 0) || (specularityCoefficient_ > 1))
    {
	FatalErrorIn
        (
            "particleThetaJohnsonJacksonFvPatchScalarField::"
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

    const fvPatchVectorField& Ua = 
	patch().lookupPatchField<volVectorField, vector>("Ua");

    const fvPatchScalarField& alpha = 
	patch().lookupPatchField<volScalarField, scalar>("alpha");

    const fvPatchScalarField& g0 = 
	patch().lookupPatchField<volScalarField, scalar>("gs0");

    const fvPatchScalarField& kappa = 
	patch().lookupPatchField<volScalarField, scalar>("kappa");
      
    scalarField magUaPatch = mag(Ua.patchInternalField());

    scalarField ThetaPatch = max(patchInternalField(), 1.0e-6);
    
    scalarField alphaPatch = max(alpha.patchInternalField(), 1.0e-6);
    
    scalarField g0Patch = g0.patchInternalField();
    
    scalarField kappaPatch = max(kappa.patchInternalField(), 1.0e-15);
 
    // The mixed BC in OpenFOAM is implemented as
    //
    // valueFraction*(Theta - ThetaRef) + (1-valueFraction)*(grad(Theta) - grad(Theta)_Ref) = 0
    //
    // To find valueFraction, we re-write Johnson and Jackson BC as
    //
    // c*(Theta - ThetaRef) + (grad(Theta) - grad(Theta)_Ref = 0
    //
    // where
    //
    // c = valueFraction/(1 - valueFraction)
    //
    // As a consequence
    //
    // valueFraction = 1/(c+1)
    //
    // We distinguish two cases, according to the value of the restituition
    // coefficient (See below)
    
    //scalarField delta = 1.0/patch().deltaCoeffs();

    if (restitutionCoefficient_ != 1.0)
    {
      	// If the restitution coefficient is < 1, Johnson and Jackson BC
	// can be written in the form
	//
	// c*(Theta - ThetaRef) + grad(Theta) = 0
	//
	// with
	//
	// c = valueFraction/(1 - valueFraction)
	//
	// grad(Theta)_Ref = 0
	//
	// and
	//
	// valueFraction = 1/(c + 1)

	this->refValue() = 2.0*specularityCoefficient_*sqr(magUaPatch)
	    /(3.0*(scalar(1) - sqr(restitutionCoefficient_)));

	this->refGrad() = 0.0;

	scalarField c = -M_PI*alphaPatch*rhoa.value()*g0Patch
	    *(scalar(1) - sqr(restitutionCoefficient_))*sqrt(3.0*ThetaPatch)
	    /(4.0*kappaPatch*alphaMax.value());

	this->valueFraction() = c/(c + scalar(1));
    }
    else
    {
	// If the restitution coefficient is 1, the BC degenerates in the form
	//
	// grad(Theta) - grad(Theta)_Ref = 0
	//
	// with 
	//
	// ThetaRef = 0
	//
	// and
	//
	// valueFraction = 0
	
	this->refValue() = 0.0;

	this->refGrad() = M_PI*specularityCoefficient_*alphaPatch*rhoa.value()
	    *g0Patch*sqrt(3.0*ThetaPatch)*sqr(magUaPatch)/(6.0*kappaPatch);

	this->valueFraction() = 0.0;
    } 

    mixedFvPatchScalarField::updateCoeffs();
}


// Write
void particleThetaJohnsonJacksonFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("restitutionCoefficient")
        << restitutionCoefficient_ << token::END_STATEMENT << nl;
    os.writeKeyword("specularityCoefficient")
        << specularityCoefficient_ << token::END_STATEMENT << nl;
    refValue().writeEntry("refValue", os);
    refGrad().writeEntry("refGradient", os);
    valueFraction().writeEntry("valueFraction", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    particleThetaJohnsonJacksonFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
