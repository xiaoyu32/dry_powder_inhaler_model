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

\*---------------------------------------------------------------------------*/

#include "meanChargeFluxFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meanChargeFluxFvPatchVectorField::
meanChargeFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    EfName_("Ef"),
    EfBool_("EfBool"),
    ThetaName_("Theta"),
    phaseName_("alpha1"),
    chargeName_("charge"),
    deltaWorkFunction_(0.0),
    ew_(1.0)     
{}


Foam::meanChargeFluxFvPatchVectorField::
meanChargeFluxFvPatchVectorField
(
    const meanChargeFluxFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    EfName_(ptf.EfName_),
    EfBool_(ptf.EfBool_),    
    ThetaName_(ptf.ThetaName_),
    phaseName_(ptf.phaseName_),
    chargeName_(ptf.chargeName_),
    deltaWorkFunction_(ptf.deltaWorkFunction_),
    ew_(ptf.ew_)  
{}


Foam::meanChargeFluxFvPatchVectorField::
meanChargeFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    EfName_(dict.lookupOrDefault<word>("Ef", "Ef")),
    EfBool_(dict.lookupOrDefault<Switch>("Efcontribution", true)),    
    ThetaName_(dict.lookupOrDefault<word>("Theta", "Theta")),
    phaseName_(dict.lookupOrDefault<word>("alpha1", "alpha1")),
    chargeName_(dict.lookupOrDefault<word>("charge", "charge")),
    deltaWorkFunction_(readScalar(dict.lookup("deltaWorkFunction"))),
    ew_(readScalar(dict.lookup("ew")))           
{
    if(!EfBool_) 
      Info << "Electric field contribution is not active on the charge transfer at the wall \n" << endl;	
}


Foam::meanChargeFluxFvPatchVectorField::
meanChargeFluxFvPatchVectorField
(
    const meanChargeFluxFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    EfName_(pivpvf.EfName_),
    EfBool_(pivpvf.EfName_),
    ThetaName_(pivpvf.ThetaName_),
    phaseName_(pivpvf.phaseName_),
    chargeName_(pivpvf.chargeName_),
    deltaWorkFunction_(pivpvf.deltaWorkFunction_),
    ew_(pivpvf.deltaWorkFunction_)
{}


Foam::meanChargeFluxFvPatchVectorField::
meanChargeFluxFvPatchVectorField
(
    const meanChargeFluxFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    EfName_(pivpvf.EfName_),
    EfBool_(pivpvf.EfName_),
    ThetaName_(pivpvf.ThetaName_),
    phaseName_(pivpvf.phaseName_),
    chargeName_(pivpvf.chargeName_),
    deltaWorkFunction_(pivpvf.deltaWorkFunction_),
    ew_(pivpvf.ew_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meanChargeFluxFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const scalar Pi = constant::mathematical::pi;
    
    const label patchI = patch().index(); 

    // cache the patch face-normal vectors
    //tmp<vectorField> nf(patch().nf());
    const vectorField nf(this->patch().nf());

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
    if(!EfBool_) Efw = vector(0,0,0); 

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
    
    //scalarField alphaMinus = alphaw * ew_ / (1.+ew_); // Numerical treatment: in the charge transport equation, div(alpha*sigmaq*Ef) 
    scalarField alphaMinus(alphaw.size(),scalar(ew_ / (1.+ew_)));
    scalarField sigmaq =   Foam::pow(2.,19./10.) * sqrtPi 
    			 * epsilon0 * alphaMinus * rhop			
			 * rstar * Foam::pow(15.0*massp/(16.0*Ystar*sqrt(rstar)),2./5.)
	      	         * 1.24127	// Gamma(12/5) 
	                 * Foam::pow(Thetaw,9./10.);
			 
    Info << "chargeBc: min(sigmaq) = " << min(sigmaq) << " max(sigmaq) = " << max(sigmaq) << endl; 			  
    
    // retrieve non-const access to charge from the database    
    volScalarField& charge =
        const_cast<volScalarField&>
        (
            db().lookupObject<volScalarField>(chargeName_)
        );
	    
    const scalar ew1 = 1./(1.+ew_);
    scalarField Efwnf(Acp.size(),scalar(0));
    Efwnf = ( Efw & nf ); 
    Info << "chargeBc: min(Efwnf) = " << min(Efwnf) << " max(Efwnf) = " << max(Efwnf) << endl;         
    scalarField qpMinus(Acp.size(),scalar(0.));  
    
    vectorField theta1 = sigmaq * Efw;
    Info << "chargeBc: min(theta1) = " << min(theta1) << " max(theta1) = " << max(theta1) << endl; 
    
    vectorField theta2 = sigmaq * deltaWorkFunction_ / ( deltaC*electronCharge ) * nf;
    Info << "chargeBc: min(theta2) = " << min(theta2) << " max(theta2) = " << max(theta2) << endl; 
            
    forAll(Acp,II)
    {
    	qpMinus[II] =   ew1 * epsilon0 * Acp[II]
                      * ( deltaWorkFunction_/(deltaC*electronCharge) + Efwnf[II] )
		      / ( 1 - ew1*max(Acp[II],1e-64)*2./(Pi*sqr(dp)) );
    }
    vectorField theta3 = - 2. * sigmaq /( Pi*epsilon0*sqr(dp) ) * qpMinus * nf;	      
    Info << "chargeBc: min(theta3) = " << min(theta3) << " max(theta3) = " << max(theta3) << endl; 
         
    vectorField::operator=(  (theta1 + theta2 + theta3)  );     
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::meanChargeFluxFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("deltaWorkFunction") << deltaWorkFunction_ << token::END_STATEMENT << nl;
    os.writeKeyword("ew") << ew_ << token::END_STATEMENT << nl;
    os.writeKeyword("Efcontribution") << EfBool_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        meanChargeFluxFvPatchVectorField
    );
}

// ************************************************************************* //
