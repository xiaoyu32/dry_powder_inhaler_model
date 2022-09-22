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

#include "meanChargeFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meanChargeFluxFvPatchScalarField::
meanChargeFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    EfName_("Ef"),
    ThetaName_("Theta"),
    phaseName_("alpha1"),
    chargeName_("charge"),
    deltaWorkFunction_(0.0),
    ew_(1.0)        
{}


Foam::meanChargeFluxFvPatchScalarField::
meanChargeFluxFvPatchScalarField
(
    const meanChargeFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    EfName_(ptf.EfName_),
    ThetaName_(ptf.ThetaName_),
    phaseName_(ptf.phaseName_),
    chargeName_(ptf.chargeName_),
    deltaWorkFunction_(ptf.deltaWorkFunction_),
    ew_(ptf.ew_)        
{}


Foam::meanChargeFluxFvPatchScalarField::
meanChargeFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    EfName_(dict.lookupOrDefault<word>("Ef", "Ef")),
    ThetaName_(dict.lookupOrDefault<word>("Theta", "Theta")),
    phaseName_(dict.lookupOrDefault<word>("alpha1", "alpha1")),
    chargeName_(dict.lookupOrDefault<word>("charge", "charge")),
    deltaWorkFunction_(readScalar(dict.lookup("deltaWorkFunction"))),
    ew_(readScalar(dict.lookup("deltaWorkFunction")))             
{}


Foam::meanChargeFluxFvPatchScalarField::
meanChargeFluxFvPatchScalarField
(
    const meanChargeFluxFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    EfName_(tppsf.EfName_),
    ThetaName_(tppsf.ThetaName_),
    phaseName_(tppsf.phaseName_),
    chargeName_(tppsf.chargeName_),
    deltaWorkFunction_(tppsf.deltaWorkFunction_),
    ew_(tppsf.deltaWorkFunction_)
{}


Foam::meanChargeFluxFvPatchScalarField::
meanChargeFluxFvPatchScalarField
(
    const meanChargeFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    EfName_(tppsf.EfName_),
    ThetaName_(tppsf.ThetaName_),
    phaseName_(tppsf.phaseName_),
    chargeName_(tppsf.chargeName_),
    deltaWorkFunction_(tppsf.deltaWorkFunction_),
    ew_(tppsf.ew_)            
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meanChargeFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }	

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const scalar Pi = constant::mathematical::pi;
    
    const label patchI = patch().index(); 

    // cache the patch face-normal vectors
    tmp<vectorField> nf(patch().nf());

    // retrieve non-const access to alpha field from the database
    volScalarField& alpha =
        const_cast<volScalarField&>
        (
            db().lookupObject<volScalarField>(phaseName_)
        );
    scalarField& alphaw = alpha.boundaryField()[patchI];
       
    // retrieve non-const access to Theta field from the database
    volScalarField& Theta =
        const_cast<volScalarField&>
        (
            db().lookupObject<volScalarField>(ThetaName_)
        );
    scalarField& Thetaw = Theta.boundaryField()[patchI];

    // retrieve non-const access to electric field from the database
    volVectorField& Ef =
        const_cast<volVectorField&>
        (
            db().lookupObject<volVectorField>(EfName_)
        );
    vectorField& Efw = Ef.boundaryField()[patchI];    

    const dictionary& electroProps =
        db().lookupObject<IOdictionary>("electrostaticsKTProperties");
	
    const scalar epsilon0 = dimensionedScalar(electroProps.lookup("epsilon0")).value(); 
    const scalar YoungModulus = dimensionedScalar(electroProps.lookup("YoungModulus")).value();            
    const scalar PoissonRatio = dimensionedScalar(electroProps.lookup("PoissonRatio")).value();
    const scalar deltaC = dimensionedScalar(electroProps.lookup("deltaC")).value();
    const scalar electronCharge = dimensionedScalar(electroProps.lookup("electronCharge")).value();

    const dictionary& transProps = db().lookupObject<IOdictionary>("transportProperties");
    const dictionary& transPhase1Props = transProps.subDict("phase1");
    const scalar dp = dimensionedScalar(transPhase1Props.lookup("d")).value(); 
    const scalar rhop = dimensionedScalar(transPhase1Props.lookup("rho")).value(); 
                
    const scalar rstar = dp/4.;
    const scalar Ystar = YoungModulus/(2.*(1-sqr(PoissonRatio)));
    const scalar massp = Pi/6.*dp*dp*dp*rhop;

    scalarField Acp =   Foam::pow(2.,12./5.) * sqrtPi 
              	      * rstar * Foam::pow(15.0*massp/(16.0*Ystar*sqrt(rstar)),2./5.)
	      	      * 0.961766	// Gamma(19/10) 
	              * Foam::pow(Thetaw,2./5.);
    
    scalarField alphaMinus = alphaw * ew_ / (1.+ew_);
    scalarField sigmaq =   Foam::pow(2.,12./5.) * sqrtPi 
    			 * epsilon0 * alphaMinus * massp
			 * rstar * Foam::pow(15.0*massp/(16.0*Ystar*sqrt(rstar)),2./5.)
	      	         * 1.24127	// Gamma(12/5) 
	                 * Foam::pow(Thetaw,9./10.); 
    
    // retrieve non-const access to charge from the database
    volScalarField& charge =
        const_cast<volScalarField&>
        (
            db().lookupObject<volScalarField>(chargeName_)
        );
    
    const scalar ew1 = 1./(1.+ew_);
    scalarField Efwnf = Efw & nf;
    scalarField dsN = Ef.mesh().Sf() & nf;  
     
    scalarField qpMinus(0,Acp.size());
    
    scalarField theta1 = sigmaq * dsN;
    scalarField theta2 = sigmaq * deltaWorkFunction_ / ( deltaC*electronCharge ) * dsN;
    scalarField theta3(0,Acp.size());
            
    forAll(Acp,II)
    {
    	qpMinus[II] =   ew1 * epsilon0 * Acp[II]
                      * ( deltaWorkFunction_/(deltaC*electronCharge) + Efwnf[II] )
		      / ( 1 - ew1*Acp[II]*2./(Pi*sqr(dp)) );
		      
	theta3[II] = 2. * sigmaq[II] /( Pi*epsilon0*sqr(dp)) * qpMinus[II] * dsN[II];	      
    }

   

        
    gradient() = ( theta1 + theta2 + theta3 );

    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        meanChargeFluxFvPatchScalarField
    );
}

// ************************************************************************* //
