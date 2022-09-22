/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "solidStresses.H"
#include "addToRunTimeSelectionTable.H"

//#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidStresses, 0);

addToRunTimeSelectionTable
(
    forceModel,
    solidStresses,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidStresses::solidStresses
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    solidVelocityFieldName_(propsDict_.lookup("solidVelocityFieldName")),
    Us_(sm.mesh().lookupObject<volVectorField> (solidVelocityFieldName_)),    
    phic_(readScalar(propsDict_.lookup("phic"))),
    interpolation_(false),
    Pp_
    (
        IOobject
        (
            "Pp",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0), 0.0)
    ),
    gradPp_
    (
        IOobject
        (
            "gradPp",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("0", dimensionSet(1,-2,-2,0,0), vector::zero)
    ),
    onlyPp_(false),
    blendedModel_(false),
    mup_
    (
        IOobject
        (
            "mup",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("0", dimensionSet(1,-1,-1,0,0), 0.0)
    ),
    gradTaup_
    (
        IOobject
        (
            "gradTaup",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("0", dimensionSet(1,-2,-2,0,0), vector::zero)
    ),
    verbose_(false),
    rhoParticle_(readScalar(propsDict_.lookup("rhoParticle"))),
    dp_(readScalar(propsDict_.lookup("dp"))),
    phid_(0.01)
{
    if (propsDict_.found("verbose")) verbose_=true;
    treatDEM_=true;
    if (propsDict_.found("interpolation"))
    {
        Info << "Solid stresses will be interpolated at the particle's position" << endl;
        interpolation_=true;
    }
    if (propsDict_.found("onlyPp"))
    {
        onlyPp_=true;
        Info << "Only solid pressure is taken into acccount" << endl;
    }

    if (propsDict_.found("blendedModel"))
    {
        blendedModel_=true;
        Info << "Using Blended Model!!!!!!!!!!!!!!!!!!" << endl;
    }

    if (propsDict_.found("alphaResidual"))
    {
        phid_ = readScalar(propsDict_.lookup("alphaResidual"));
        Info << "Using Blended Model!!!!!!!!!!!!!!!!!!" << endl;
    }    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

solidStresses::~solidStresses()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidStresses::setForce() const
{
    scalar alpp(0.);
    label cellI(-1);
    vector position(0,0,0);
    scalar voidfraction(1);
    vector gradPp(0,0,0);
    vector gradTaup(0,0,0);    
    scalar ds(0);
    scalar Vp(0);
    vector force(0,0,0);   

    volTensorField gradUs = fvc::grad(Us_);
    volTensorField gradUsT = fvc::grad(Us_)().T();
    volTensorField strainRate = gradUs + gradUsT - ((2.0/3.0)*I)*tr(gradUs);
    volScalarField SS = strainRate && strainRate;
    
    scalar constB1  =  0.628;
    scalar constC1  = -1.073;    
    scalar constChi =  1.650;

    scalar constB2   =  0.351;
    scalar constC2   = -0.895;
    scalar constChi2 =  1.650;   
        
    //scalar phid = 0.01;
    scalar phia = phic_ - phid_;
    scalar k = 700.0;
    scalar small = 1e-12;
    scalar deltaAlpp(0.);
    forAll(Pp_,cellI)
    {
        alpp = 1.-voidfraction_[cellI];
	
	// Alpha diff. from phia
	deltaAlpp = alpp - phia; 
	
		    if(blendedModel_)
	 	    {


			    if( alpp < phia )
			    {
			       Pp_[cellI]  =                             pow(
			    pow(   0.015 * SS[cellI]*dp_*dp_*rhoParticle_/pow(max(small,sqrt((alpp-phic_)* (alpp-phic_))),2.14)  ,-1.0)
			    +
                            pow(k* 0.1/dp_ * pow(sqrt(SS[cellI])*dp_/sqrt(k/(rhoParticle_*dp_)),0.477),-1.0)
			    , -1.0);
			    }else
			    {

			       	
			       Pp_[cellI]  =  k*0.2/dp_ * pow(max(small,sqrt((alpp-phic_)* (alpp-phic_))),0.67) + k* 0.1/dp_ * pow(sqrt(SS[cellI])*dp_/sqrt(k/ (rhoParticle_*dp_)),0.477);
			    }
		    }else
		    {
		    
			    if( alpp < phia )
			    {
			       Pp_[cellI]  = ( constB1*pow(alpp,2.50)+constC1*pow(alpp,3.5) )
	        	        	    /( pow(phic_-alpp,constChi) 		    ) 
			        	    * rhoParticle_ * dp_ * dp_ 
			        	    * SS[cellI];
			    }else
			    {

			       	
			       Pp_[cellI]  = (  ( constB1*pow(phia,2.50)+constC1*pow(phia,3.5) )
	        		               /( pow(phic_-phia,constChi ) 		       )
					      + deltaAlpp
					       *(
						 (1.65 * (constB1*pow(phia,2.5)+      constC1*pow(phia,3.5)))/pow(phic_-phia,2.65)+
						 (2.50 *  constB1*pow(phia,1.5)+3.5 * constC1*pow(phia,2.5) )/pow(phic_-phia,1.65)
						 )
					      )* rhoParticle_ * dp_ * dp_ 
				               * SS[cellI];
			    }		    
		    } 	

		    if( alpp < phia )
		    { 
		       mup_[cellI] = ( constB2*pow(alpp,1.75)+constC2*pow(alpp,3.5) )
	                            /( pow(phic_-alpp,constChi2) )
		                    * rhoParticle_ * dp_ * dp_ 
		                    *sqrt(SS[cellI]);			    		    
		    }else
		    {
		       mup_[cellI] = (  ( constB2*pow(phia,1.75)+constC2*pow(phia,3.5) )
	                               /( pow(phic_-phia,constChi2) 		       )
				      + deltaAlpp
				       *(
				         (1.65 * (constB2*pow(phia,1.75)+       constC2*pow(phia,3.5)))/pow(phic_-phia,2.65)+
				         (1.75 *  constB2*pow(phia,0.75)+ 3.5 * constC2*pow(phia,2.5) )/pow(phic_-phia,1.65)
				        )
				      )* rhoParticle_ * dp_ * dp_ 
				       *sqrt(SS[cellI]);			   		    
		    }


 //	alpp = min(alpp,phic_-0.001);
 //       Pp_[cellI]  = ( constB1*pow(alpp,2.50)+constC1*pow(alpp,3.5) )
//	             /( pow(phic_-alpp,constChi) )
//		     * rhoParticle_ * dp_ * dp_ 
//		     * SS[cellI]+10000.0/dp_*pow(max(0.0,alpp-phic_),1.10);
//	mup_[cellI] = ( constB2*pow(alpp,1.75)+constC2*pow(alpp,3.5) )
//	             /( pow(phic_-alpp,constChi2) )
//		     * rhoParticle_ * dp_ * dp_ 
//		     *sqrt(SS[cellI]);		     	      
    }

    gradPp_ = -fvc::grad(Pp_);   
    gradTaup_ = fvc::div(mup_*strainRate);
    
    interpolationCellPoint<vector> gradPpInterpolator_(gradPp_);
    interpolationCellPoint<vector> gradgradTaupInterpolator_(gradTaup_);    
    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            force = vector(0,0,0);

            if (cellI > -1) // particle Found
            {
                if(interpolation_)
                {
	            position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    gradPp = gradPpInterpolator_.interpolate(position,cellI);
		    gradTaup = gradgradTaupInterpolator_.interpolate(position,cellI);
                }else
                {
                    voidfraction = voidfraction_[cellI];
                    gradPp = gradPp_[cellI];
		    gradTaup = gradTaup_[cellI];
                }
                ds = 2*particleCloud_.radius(index);
                Vp = ds*ds*ds*M_PI/6;

                // calc particle's force
                if (voidfraction<1)
		{
                    if(onlyPp_)
		    {
		      force = gradPp/(1-voidfraction)*Vp;
		    }else
		    {
		      force = (gradPp+gradTaup)/(1-voidfraction)*Vp;
		    } 	
		}
            }

            if(verbose_ && index ==0 )
            {
                Info << "cellI = " << cellI << endl;
                Info << "voidfraction = " << voidfraction << endl;
                Info << "Grad. of solid pressure = " << gradPp << endl;
		Info << "Div. of solid shear = " << gradTaup << endl;
            }

            // set force on particle
            if(!treatDEM_){
                if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += force[j];
                else  for(int j=0;j<3;j++) impForces()[index][j] += force[j];
            }
            for(int j=0;j<3;j++) DEMForces()[index][j] += force[j];
        //}
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
