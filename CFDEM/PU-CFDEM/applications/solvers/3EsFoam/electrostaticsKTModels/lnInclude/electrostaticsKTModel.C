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

#include "electrostaticsKTModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electrostaticsKTModel::electrostaticsKTModel
(
    const phaseModel& phase1,
    const volScalarField& alpha1
)
:
    phase1_(phase1),
    U1_(phase1.U()),
    alpha1_(alpha1),
    rho1_(phase1.rho()),
    da_(phase1.d()),
    electrostaticsKTProperties_
    (
        IOobject
        (
            "electrostaticsKTProperties",
            U1_.time().constant(),
            U1_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    electrostaticsKT_(electrostaticsKTProperties_.lookup("electrostaticsKT")),
    epsilon0_(electrostaticsKTProperties_.lookup("epsilon0")),
    YoungModulus_(electrostaticsKTProperties_.lookup("YoungModulus")),            
    PoissonRatio_(electrostaticsKTProperties_.lookup("PoissonRatio")),            
    kineticTheoryProperties_
    (
        IOobject
        (
            "kineticTheoryProperties",
            U1_.time().constant(),
            U1_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    alphaMax_(kineticTheoryProperties_.lookup("alphaMax")),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            kineticTheoryProperties_
        )
    ),
    gs0E_
    (
        IOobject
        (
            "gs0E",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    Theta_
    (
        IOobject
        (
            "Theta",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh().lookupObject<volScalarField> ("Theta")
    ),        
    sigmaq_
    (
        IOobject
        (
            "sigmaq",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, -3, 3, 0, 0, 2, 0), 0.0)
    ),
    kappaq_
    (
        IOobject
        (
            "kappaq",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
    ),
    chargeTurbulentDiffusivity_(electrostaticsKTProperties_.lookup("chargeTurbulentDiffusivity")),
    kappaqt_
    (
        IOobject
        (
            "kappaqt",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
    ),
    mupt_(dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)),
    DTheta_
    (
        IOobject
        (
            "DTheta",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
    ),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    e_(kineticTheoryProperties_.lookup("e")),
    mu1E_
    (
        IOobject
        (
            "mu1E",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
    ),
    sigmaqEf_
    (
        IOobject
        (
            "sigmaqEf_",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedVector("zero",  dimensionSet(1, -2, 0, 0, 0, 1, 0), vector(0.0,0.0,0.0))
    )                       
{
    Info << "Electrostatics KT is " << electrostaticsKT_ << endl;
    Info << "Charge diffusivity by phase fluctuations is " << chargeTurbulentDiffusivity_ << endl;
    if(chargeTurbulentDiffusivity_) 
    {
    	mupt_ = dimensionedScalar(electrostaticsKTProperties_.lookup("solidPhaseTurbulentViscosity"));		
    }	
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electrostaticsKTModel::compute()
{
    
    if (!electrostaticsKT_)
    {
        return;
    }

    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    // 1/rstar = 1/r1 + 1/r2
    dimensionedScalar rstar("zero",dimensionSet(0,1,0,0,0,0,0),da_.value()/4.);
    
    
    // 1/Y = (1-v1^2)/Y1 + (1-v2^2)/Y2 
    dimensionedScalar Ystar
    (
    	"zero",
	dimensionSet(1,-1,-2,0,0,0,0),
	YoungModulus_.value()/(2.*(1-sqr(PoissonRatio_.value())))
    );  

    // particle mass
    dimensionedScalar massp
    (
    	"zero",
	dimensionSet(1, 0, 0, 0, 0, 0, 0),
	constant::mathematical::pi/6.*da_.value()*da_.value()*da_.value()*rho1_.value()
    );
    
    // Calculating the radial distribution function 
    gs0E_ = radialModel_->g0
    (
        min(max(alpha1_, scalar(1e-6)), alphaMax_ - 0.01),
        alphaMax_
    );

    // Calculating tribocharging conductivity (m_p*d_p^3*n_p^2 = 1./6.*pi*alpha_p^2*rho_p)
    sigmaq_ =   pow(2.,9./5.) * sqrtPi * epsilon0_ * gs0E_
              //* sqr(alpha1_) * rho1_ 
	      * alpha1_ * rho1_ 	// Numerical treatment: in the charge transport equation, div(alpha*sigmaq*Ef) 
	      * 1.242169	// Gamma(12/5) 
	      * rstar 
	      * pow(15.0*massp/(16.0*Ystar*sqrt(rstar)),0.4)
	      * pow(Theta_,0.9);

    // Calculating tribocharging diffusivity 
    kappaq_ =   pow(2.,-1./5.) / sqrtPi * gs0E_
              //* sqr(alpha1_) * rho1_ / da_
	      * alpha1_ * rho1_ / da_	// Numerical treatment: in the charge transport equation, div(alpha*kappaq*grad(Qp) 
	      * 1.242169	// Gamma(12/5) 
	      * rstar 
	      * pow(15.0*massp/(16.0*Ystar*sqrt(rstar)),0.4)
	      * pow(Theta_,0.9);

    // retrieve non-const access to electric field from the database
    volVectorField& Ef =
        const_cast<volVectorField&>
        (
            U1_.mesh().lookupObject<volVectorField>("Ef")
        );
		      
    sigmaqEf_ = sigmaq_ * Ef;	      
          
    Info<< "Electrostatics KT: min(Theta) = " << min(Theta_).value()
        << ", max(Theta) = " << max(Theta_).value() << endl;
    Info<< "Electrostatics KT: min(sigmaq) = " << min(sigmaq_).value()
        << ", max(sigmaq) = " << max(sigmaq_).value() << endl;
    Info<< "Electrostatics KT: min(sigmaqEf) = " << min(sigmaqEf_).value()
        << ", max(sigmaqEf) = " << max(sigmaqEf_).value() << endl;
    Info<< "Electrostatics KT: min(kappaq) = " << min(kappaq_).value()
        << ", max(kappaq) = " << max(kappaq_).value() << endl;

    if(chargeTurbulentDiffusivity_) 
    {
	// Calculating tribocharging turbulent mass diffusivity
	DTheta_ = rho1_ * 1./(3.*sqrtPi) * sqrt(Theta_) * da_ / max(alpha1_, scalar(1e-6));
	Info<< "Electrostatics KT: min(Diff) = " << min(DTheta_).value()
            << ", max(Diff) = " << max(DTheta_).value() << endl;
	    
   	// Calculating solid viscosity
	mu1E_ = viscosityModel_->mu1(alpha1_, Theta_, gs0E_, rho1_, da_, e_);
	Info<< "Electrostatics KT: min(mu1) = " << min(mu1E_).value()
            << ", max(mu1) = " << max(mu1E_).value() << endl;
	  
	// Calculating tribocharging turbulent diffusivity    	
	dimensionedScalar mu1Esmall
	(
            "small",
            dimensionSet(1 , -1 ,-1 ,0 , 0, 0, 0),
            1.0e-6
	);
	kappaqt_ = mupt_ / ( mu1E_ + mu1Esmall ) * DTheta_ ; // (nu_p^t * rhop * D / mu_p)

	Info<< "Electrostatics KT: min(kappaqt) = " << min(kappaqt_).value()
            << ", max(kappaqt) = " << max(kappaqt_).value() << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electrostaticsKTModel::~electrostaticsKTModel()
{}

// ************************************************************************* //

