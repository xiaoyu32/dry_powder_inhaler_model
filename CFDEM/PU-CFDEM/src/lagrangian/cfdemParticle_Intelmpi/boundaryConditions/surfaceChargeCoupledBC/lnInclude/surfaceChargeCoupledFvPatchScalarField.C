/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "surfaceChargeCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from patch and internal field

surfaceChargeCoupledFvPatchScalarField::
surfaceChargeCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    nbrFieldName_("phie"),
    surfaceChargeFieldName_("surfaceChargeDensity"),
    epsilon_
    (
        dimensionedScalar("0", dimensionSet(0,2,0,0,0,0,0), 8.85419e-12 )
    ),
	epsDeltaMapped_(false)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}

//- Construct by mapping given
//  surfaceChargeCoupledFvPatchScalarField onto a
//  new patch

surfaceChargeCoupledFvPatchScalarField::
surfaceChargeCoupledFvPatchScalarField
(
    const surfaceChargeCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    nbrFieldName_(ptf.nbrFieldName_),
    surfaceChargeFieldName_(ptf.surfaceChargeFieldName_),
    epsilon_(ptf.epsilon_),
    epsDeltaMapped_(false)
{}

//- Construct from patch, internal field and dictionary
surfaceChargeCoupledFvPatchScalarField::
surfaceChargeCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    nbrFieldName_(dict.lookup("nbrField")),
    surfaceChargeFieldName_(dict.lookup("surfaceChargeField")),
    epsDeltaMapped_(false),
    epsilon_
    (
        dimensionedScalar("0", dimensionSet(0,2,0,0,0,0,0), 8.85419e-12 )
    )
{

	const fvMesh& mesh = patch().boundaryMesh().mesh();

	// Read permittivity from dictionary
	if (mesh.foundObject<IOdictionary>("physicalProperties"))
	{

		epsilon_=
			mesh.lookupObject<IOdictionary>
			(
			   "physicalProperties"
			).lookupOrDefault("epsilon0", dimensionedScalar("0", dimensionSet(0,2,0,0,0,0,0), 8.85419e-12 ) );
	} 
	else
	{
		IOdictionary electrostaticProperties
		(
			IOobject
			(
				"physicalProperties",
				mesh.time().constant(),
				mesh,
				IOobject::MUST_READ_IF_MODIFIED,
				IOobject::NO_WRITE
			)
		);

		epsilon_ = electrostaticProperties.lookupOrDefault( "epsilon0", dimensionedScalar("0", dimensionSet(0,2,0,0,0,0,0), 8.85419e-12 ) );
	}



    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
	
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    epsDelta_ = epsilon()*patch().deltaCoeffs();


    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}

//- Construct as a copy

surfaceChargeCoupledFvPatchScalarField::
surfaceChargeCoupledFvPatchScalarField
(
    const surfaceChargeCoupledFvPatchScalarField& SCCSF,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(SCCSF, iF),
    nbrFieldName_(SCCSF.nbrFieldName_),
    surfaceChargeFieldName_(SCCSF.surfaceChargeFieldName_),
    epsilon_(SCCSF.epsilon_),
	epsDeltaMapped_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void surfaceChargeCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>(patch().patch());
    //const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());

    //const label samplePatchi = mpp.samplePolyPatch().index();
    //const fvPatch& nbrPatch = refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

	// Obtain the field of the neighbouring mesh
    //const volScalarField& nbrField = 
	//	nbrPatch.lookupPatchField<volScalarField>(nbrFieldName_);
    const volScalarField& nbrField = nbrMesh.lookupObject<volScalarField>(nbrFieldName_);

    //  Calculate the electric potential at the boundary
	scalar boundaryFieldSize = nbrField.boundaryField().size();
    scalarField nbrMappedIntFld(boundaryFieldSize);

    //tmp<scalarField> nbrMappedIntFld(new scalarField(boundaryFieldSize, 0.0));

	scalarField surfCharge =
	    patch().lookupPatchField<volScalarField, scalar>
	    (
	        surfaceChargeFieldName_
	    );

	const volScalarField& nbrSurfChargeField =
		nbrMesh.lookupObject<volScalarField>
	    (
	        surfaceChargeFieldName_
	    );


    switch (mpp.mode())
    {
        default:
        {
	    

            break;
        }
        case mappedPatchBase::NEARESTPATCHFACE:
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mpp.samplePatch());

            const surfaceChargeCoupledFvPatchScalarField& nbrPatchField = 
				refCast<const surfaceChargeCoupledFvPatchScalarField>
				(nbrField.boundaryField()[nbrPatchID]);
            //const scalarField& nbrPatchInternalField = nbrPatchField.patchInternalField();

    		nbrMappedIntFld = nbrPatchField.patchInternalField();
            mpp.distribute(nbrMappedIntFld);

			if(!epsDeltaMapped_)
			{
		   		nbrEpsDelta_ = nbrPatchField.epsDelta();
    			mpp.distribute(nbrEpsDelta_);
			}

			scalarField nbrSurfCharge = nbrSurfChargeField.boundaryField()[nbrPatchID];
				
			mpp.distribute(nbrSurfCharge);
			surfCharge += nbrSurfCharge;

            break;
        }
        case mappedPatchBase::NEARESTFACE:
        {
			// Assign all the nbr boundary field values, leaving the internal values zero
            scalarField allNbrIntValues(nbrMesh.nFaces(), 0.0);
            scalarField allNbrSurfCharge(nbrMesh.nFaces(), 0.0);

            forAll(nbrField.boundaryField(), patchi)
            {
		        const fvPatchScalarField& nbrPatchField = (nbrField.boundaryField()[patchi]);
                const scalarField& nbrPatchInternalField = nbrPatchField.patchInternalField();

				const scalarField& nbrSurfCharge = nbrSurfChargeField.boundaryField()[patchi];

                label faceStart = nbrPatchField.patch().start();

                forAll(nbrPatchField, facei)
                {
                    allNbrIntValues[faceStart + facei] = nbrPatchInternalField[facei];
                    allNbrSurfCharge[faceStart + facei] = nbrSurfCharge[facei];
                }
            }
            mpp.distribute(allNbrIntValues);
            nbrMappedIntFld.transfer(allNbrIntValues);

            mpp.distribute(allNbrSurfCharge);
    		scalarField nbrMappedSurfCharge;
            nbrMappedSurfCharge.transfer(allNbrSurfCharge);
			surfCharge += nbrMappedSurfCharge;

			if(!epsDeltaMapped_)
			{

				const scalar nbrEpsilon = 
				nbrMesh.lookupObject<IOdictionary>
				(
				   "physicalProperties"
				).lookupOrDefault("epsilon0",dimensionedScalar("0", dimensionSet(0,2,0,0,0,0,0), 8.85419e-12 ) ).value();

            	scalarField allNbrEpsDelta(nbrMesh.nFaces(), 0.0);
		        forAll(nbrField.boundaryField(), patchi)
		        {
				    const fvPatchScalarField& nbrPatchField = nbrField.boundaryField()[patchi];
		            const scalarField& nbrPatchEpsDelta = nbrEpsilon*nbrPatchField.patch().deltaCoeffs();

		            label faceStart = nbrPatchField.patch().start();

		            forAll(nbrPatchField, facei)
		            {
		                allNbrEpsDelta[faceStart + facei] = nbrPatchEpsDelta[facei];
		            }
		        }
    			mpp.distribute(allNbrEpsDelta);
           	 	nbrEpsDelta_.transfer(allNbrEpsDelta);
			}

		}
	}




	/*

     When passing a surface with an surface charge of sigma, there is a voltage jump at the 
	 surface such that
	
	 	(Sf/magSf)**(myEps*myE - nbrEps*nbrE)  =  sigma,
	
	 where myE and nbrE are the electric field vectors nearby the passed surface 
    (at the neighbouring cells) and myEps and nbrEps are the respective electric permittivities. 
	The electric field at the boundary is obtained from the negative gradient of the electric
	potential phi, so the elecric field strength in the direction of the surface normal is 

		(Sf/magSf)**myE = - (myPhi -boundaryPhi)*myDelta,

	where phi_my is the potential at the cell centroid and phi_f is the potential at the 
	boundary face. Respectively

		(Sf/magSf)**nbrE = - (boundaryPhi -nbrPhi)*nbrDelta.

	This results in

	 	-myEps*(myPhi -boundaryPhi)*myDelta + eps_nbr*(boundaryPhi -nbrPhi)*nbrDelta = sigma,
		-> boundaryPhi*(myEps*myDelta + nbrEps*nbrDelta) = myPhi*myEps*myDelta + nbrPhi*nbrEps*nbrDelta + sigma
		-> boundaryPhi  = (myPhi*myEpsDelta + nbrPhi*nbrEpsDelta + sigma) / (myEpsDelta + nbrEpsDelta)

	In the mixed boundary condition this is obtained by setting

        - refValue = nbrPhi.

        - refGrad = sigma / myEps

        - valueFraction = nbrEpsDelta / (nbrEpsDelta + myEpsDelta)

		- Now the result basing on refValue is weighted by nbrEpsDelta / (nbrEpsDelta + myEpsDelta) and 
			the result basing on refGrad is weighted by myEpsDelta / (nbrEpsDelta + myEpsDelta). 

		The fixedValue bc would give

			boundaryPhi = nbrPhi

		and the fixedGradient bc would give

			refGrad = (boundaryPhi-myPhi)*myDelta

			-> boundaryPhi 	= myPhi + refGrad / myDelta
							= myPhi + sigma / myEpsDelta

		weighed by the valueFraction these two give the desired result.

	*/

    this->refValue() =  nbrMappedIntFld;

    this->refGrad() = surfCharge / epsilon();

    this->valueFraction() = nbrEpsDelta_ / (nbrEpsDelta_ + epsDelta());

	
    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << nbrMesh.name() << ':'
            << " Potential"
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void surfaceChargeCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nbrField")<< nbrFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("surfaceChargeField")<< surfaceChargeFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("epsilon0")<< epsilon_
        << token::END_STATEMENT << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    surfaceChargeCoupledFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
