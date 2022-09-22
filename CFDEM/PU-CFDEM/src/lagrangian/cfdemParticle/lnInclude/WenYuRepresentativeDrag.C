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

#include "WenYuRepresentativeDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

//#include "mpi.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WenYuRepresentativeDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    WenYuRepresentativeDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
WenYuRepresentativeDrag::WenYuRepresentativeDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    granVelFieldName_(propsDict_.lookup("granVelFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    Us_(sm.mesh().lookupObject<volVectorField> (granVelFieldName_)),      
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),   
    verbose_(false),
    interpolation_(false),
    parcelApproach_(false),
    npart_(1),
    dpart_(1.0),
    rhopart_(1.0),
    apiType_
    (
        propsDict_.lookupOrDefault("apiType", 1)
    ),
    velslipdiff_(false),
    weighting_(false),
    interpolation_void_(false),
    interpolation_vel_(false)
{
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    if (propsDict_.found("interpolation")) 
    {
        interpolation_=true;
        Info << "Using interpolated values." << endl;
    }
    
    if (propsDict_.found("interpolation_void")) 
    {
        interpolation_void_=true;
        Info << "Using interpolated voidage values." << endl;
    }
    
    if (propsDict_.found("interpolation_vel")) 
    {
        interpolation_vel_=true;
        Info << "Using interpolated velocity values." << endl;
    }
    
    // Representative particle Approach
    parcelApproach_=true;
    npart_  = readScalar(propsDict_.lookup("npart"));
    dpart_  = readScalar(propsDict_.lookup("dpart"));
    rhopart_= readScalar(propsDict_.lookup("rhopart"));   

    Info << "Using Representative Particle Approach" << endl;

	    Info << "Number of particles in a parcel = " << npart_    << endl;
	    Info << "Particle diameter               = " << dpart_    << endl;
	    Info << "Particle density 		 = " << rhopart_  << endl;
	    Info << "Parcel Type			 : " << apiType_  << endl;			

    
    // Calculation of local relative velocity at grid location, then interpolate to particle position: int<U-Us>@x_p
    if (propsDict_.found("velslipdiff")) 
    {
       velslipdiff_=true;
       interpolation_=true; // Force calculate at the particle position
       Info << "Calculation of local relative velocity at grid location by using Eulerian velocities then interpolate to particle position" << endl;  	
    }    

    if (propsDict_.found("weighting")) 
    {
    	weighting_= true;    
	interpolation_= false;
        Info << " Weighting is using for interpolation " << endl;
    }
                 
    if (propsDict_.found("implDEM"))
    {
        treatExplicit_=false;
        implDEM_=true;
        setImpDEMdrag();
        Info << "Using implicit DEM drag formulation." << endl;
    }
    particleCloud_.checkCG(false);
//
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

WenYuRepresentativeDrag::~WenYuRepresentativeDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WenYuRepresentativeDrag::setForce() const
{
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label  cellI=0;

    vector Up(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar CD(0);
    scalar Volp(0);
    scalar alps(0);
    scalar betaP(0);

    // Particle characteristics 
    scalar v_t(0); 

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    // Slip velocity 
    interpolationCellPoint<vector> UsInterpolator_(Us_);	
    
    // To calculate domain-averaged slip velocity --> averaging particle velocity
    //vector sumUp(0,0,0);
    
    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);
            betaP = 0;
            Volp = 0;
	    
	    voidfraction = 0;
            Ufluid = vector(0,0,0);
	    
            if (cellI > -1) // particle Found
            {
                if(interpolation_)
                {
	            position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;
                    if(voidfraction>1.00) voidfraction = 1.00;
                    if(voidfraction<0.40) voidfraction = 0.40;
                }
		else if(interpolation_void_)
                {
	            position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = U_[cellI];//UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;
                    if(voidfraction>1.00) voidfraction = 1.00;
                    if(voidfraction<0.40) voidfraction = 0.40;
                }
		else if(interpolation_vel_)
                {
	            position = particleCloud_.position(index);
                    voidfraction = particleCloud_.voidfraction(index);//voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;
                    if(voidfraction>1.00) voidfraction = 1.00;
                    if(voidfraction<0.40) voidfraction = 0.40;
                }
		else if(weighting_)
		{	
                    scalar sumWeighting(0);   	    
		    for (int subCell = 0; subCell < particleCloud_.voidFractionM().cellsPerParticle()[index][0]; subCell++)	
		    {
		    	label subCellID = particleCloud_.cellIDs()[index][subCell]; 
		        sumWeighting += particleCloud_.weighting(index,subCell); 
                	voidfraction += particleCloud_.weighting(index,subCell)*voidfraction_[subCellID];
			Ufluid       += particleCloud_.weighting(index,subCell)*voidfraction_[subCellID]*U_[subCellID];
		    } 
		    Ufluid /= voidfraction;	
                    //voidfraction /= sumWeighting; 
                    //if(verbose_ && index==0) Info << "sumWeighting = " << sumWeighting << endl;  
		}
		else
                {
		    voidfraction = particleCloud_.voidfraction(index);
                    Ufluid = U_[cellI];
                }

                Up = particleCloud_.velocity(index);
                Ur = Ufluid-Up;
                ds = 2*particleCloud_.radius(index);

		// DomainAveVelSlip
		//if(DomainAveVelSlip_) for(int j=0;j<3;j++) sumUp[j] += Up[j];


		if(velslipdiff_)
 		{
    		    	Ur = UInterpolator_.interpolate(position,cellI)-UsInterpolator_.interpolate(position,cellI);	
		}
		
		// Representative approach: drag based on particle diameter
		int ptype = particleCloud_.types()[index][0];
		//Info <<"For partice index "<<index<<" ptype is "<<ptype<<endl;
		//Info <<"apiType_ is "<<apiType_<<endl;
		
		if ( (parcelApproach_) && (ptype == apiType_) )
		{
	 	    ds = dpart_;
		}		
                
		nuf = nufField[cellI];
                rho = rho_[cellI];
                magUr = mag(Ur);
                Volp = ds*ds*ds*M_PI/6;
                alps = 1-voidfraction+SMALL;

                if (magUr > 0)
                {
                    // calc particle Re number
		      Rep = voidfraction*ds*magUr/(nuf+SMALL);
		      
                    // calc CD
		    if (Rep < 1000)
		    {
		      CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-1.65);
		    }else
		    {
                      CD = 0.44*pow(voidfraction,-1.65); 
		    }  
		     
		    // calc drag coefficient 
		      betaP = 3./4.*rho*CD*magUr/ds;  	

                    // calc particle's drag
                    drag = Volp*betaP*Ur;
			
		    
		    // Representative particle approach: drag on one particle multiply by number of particles in a parcel
		    if ( (parcelApproach_) && (ptype == apiType_) )
		    {
			    //Info << "Using WenYuRepresentativeDrag ..."<<endl;
			    drag *= npart_;
		    }	
		      
                    if (modelType_=="B")
                        drag /= voidfraction;
                }

                if(verbose_ && index <1)
                {
                    Info << "index = " << index << endl;
                    Info << "Up = " << Up << endl;
                    Info << "Ur = " << Ur << endl;
                    Info << "dp = " << ds << endl;
                    Info << "dparcel = " << 2*particleCloud_.radius(index) << endl;
                    Info << "rho = " << rho << endl;
                    Info << "nuf = " << nuf << endl;
                    Info << "voidfraction = " << voidfraction << endl;
                    Info << "Rep = " << Rep << endl;
                    Info << "Volp = " << Volp << endl;
                    Info << "alps = " << alps << endl;
                    Info << "CD = " << CD << endl;
		    Info << "betaP = " << betaP << endl;		   		    
                    Info << "drag = " << drag << endl;
                }
            }
	   
	// Mostafa
	// removied the + sign befor drag[j] and drag[j]    	    
            // set force on particle
      //      Info << "------------------No force over addition -------------------"<<endl;
            if(treatExplicit_) 
	    {
	    	for(int j=0;j<3;j++) expForces()[index][j] += drag[j];
            }
	    else
	    {
	    	for(int j=0;j<3;j++) impForces()[index][j] += drag[j];
	    }
	    
            // set Cd
            if(implDEM_)
            {
                for(int j=0;j<3;j++) fluidVel()[index][j]=Ufluid[j];

                if (modelType_=="B")
                    Cds()[index][0] = Volp*betaP/voidfraction;
                else
                    Cds()[index][0] = Volp*betaP;	    	
    	
            }
	    else
	    {
                for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];
            }


    }
	/*
	if(DomainAveVelSlip_)
	{
		Info << "Domain-averaged slip velocity = " 
	     << fvc::domainIntegrate(voidfraction_*U_).value()/fvc::domainIntegrate(voidfraction_).value()
			  -sumUp/particleCloud_.numberOfParticles() << endl;

	}    

	if(DomainAveMom_)
	{
		Info 	<< "Domain-averaged fluid momentum = " 
	     		<< fvc::domainIntegrate(rho_*voidfraction_*U_).value()
			<< " total particle Vol*Up = " << sumUp*ds*ds*ds*M_PI/6 << endl;	
	}
	*/
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
