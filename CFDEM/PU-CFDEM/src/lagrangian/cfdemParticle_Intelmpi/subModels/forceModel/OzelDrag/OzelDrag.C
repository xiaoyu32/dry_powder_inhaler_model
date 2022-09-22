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

#include "OzelDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

//#include "mpi.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(OzelDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    OzelDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
OzelDrag::OzelDrag
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
    Kxx_(readScalar(propsDict_.lookup("Kxx"))),
    Kyy_(readScalar(propsDict_.lookup("Kyy"))),
    Kzz_(readScalar(propsDict_.lookup("Kzz"))),    
    velslipdiff_(false),     
        gravityFieldName_(propsDict_.lookup("gravityFieldName")),
    #if defined(version21) || defined(version16ext)
        g_(sm.mesh().lookupObject<uniformDimensionedVectorField> (gravityFieldName_))
    #elif defined(version15)
        g_(dimensionedVector(sm.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(gravityFieldName_)).value())
    #endif
{
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    if (propsDict_.found("interpolation")) 
    {
        interpolation_=true;
        Info << "Using interpolated value of U." << endl;
    }
    
    // Parcel Approach
    if (propsDict_.found("parcel"))
    {
    
        parcelApproach_=true;
        npart_  = readScalar(propsDict_.lookup("npart"));
        dpart_  = readScalar(propsDict_.lookup("dpart"));
        rhopart_= readScalar(propsDict_.lookup("rhopart"));   
	 
        Info << "Using Parcel Approach" << endl;

	        Info << "Number of particles in a parcel = " << npart_    << endl;
		Info << "Particle diameter               = " << dpart_    << endl;
		Info << "Particle density 		 = " << rhopart_  << endl;
			
    }
    
    // Calculation of local relative velocity at grid location, then interpolate to particle position: int<U-Us>@x_p
    if (propsDict_.found("velslipdiff")) 
    {
       velslipdiff_=true;
       interpolation_=true; // Force calculate at the particle position
       Info << "Calculation of local relative velocity at grid location by using Eulerian velocities then interpolate to particle position" << endl;  	
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

OzelDrag::~OzelDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void OzelDrag::setForce() const
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

    // Particle characteristics 
    scalar grav(0);
    scalar v_t(0); 

    // Sub-grid variables and constants
    scalar rhos(0);
    scalar taup(0);

    scalar func_f(0);
    scalar delta_bar(0);
    scalar const_f1(0.15);

    scalar func_h(0);
    scalar const_h1(0.10);
    scalar const_h2(-1.88); 
    scalar const_h3(5.16);

    const scalar alpsmax(0.64);
 
    // Kronecker delta
    tensor deltaij(1,0,0,0,1,0,0,0,1);
    tensor Kii(0,0,0,0,0,0,0,0,0);
    // Drag coefficient
    tensor betaP(0,0,0,0,0,0,0,0,0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    // Slip velocity 
    interpolationCellPoint<vector> UsInterpolator_(Us_);	
    
    // Gravity
    grav = mag(g_.value());	

   
    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
    	    // Drag correction tensor
      		// Kxx(3.3);
      		// Kyy(3.3);
      		// Kzz(5.0);
            Kii = tensor(Kxx_,0,0,0,Kyy_,0,0,0,Kzz_);
            drag = vector(0,0,0);
            betaP = tensor(0,0,0,0,0,0,0,0,0);
            Volp = 0;
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
                }else
                {
		    voidfraction = particleCloud_.voidfraction(index);
                    Ufluid = U_[cellI];
                }

                Up = particleCloud_.velocity(index);
                Ur = Ufluid-Up;

		if(velslipdiff_)
 		{
    		    	Ur = UInterpolator_.interpolate(position,cellI)-UsInterpolator_.interpolate(position,cellI);	
		}
				
                ds = 2*particleCloud_.radius(index);
		rhos = rhopart_;		
		
		// Parcel approach: drag based on particle diameter
		if(parcelApproach_)
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
		      CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-2.65);
		    }else
		    {
                      CD = 0.44*pow(voidfraction,-2.65); 
		    }   

		    // calc taup
		      taup = 4./3.*rhos/rho*ds/(CD*magUr);

		    // Mesh dependency of drag correction
		    // Dimensionless characteristic mesh length
		      delta_bar = pow(particleCloud_.mesh().V()[cellI],1./3.)/(taup*magUr);
		      func_f = pow(delta_bar,2)/(const_f1+pow(delta_bar,2));	 		
		     
                    // Volume fraction dependency of drag correction
		      func_h = - tanh(alps/const_h1)*sqrt(alps/alpsmax)*pow((1-alps/alpsmax),2)
                                *(1+const_h2*alps/alpsmax+const_h3*pow(alps/alpsmax,2));

		    // Limiting Kii*func_h*func_f <=-0.999	
			if( Kii.xx()*func_h*func_f <=-0.999 )
			{
				Kii.xx() = -0.999/(func_h*func_f);
			}
 
                        if( Kii.yy()*func_h*func_f <=-0.999 )
                        {
                                Kii.yy() = -0.999/(func_h*func_f);
                        }

                        if( Kii.zz()*func_h*func_f <=-0.999 )
                        {
                                Kii.zz() = -0.999/(func_h*func_f);
                        }

		    // calc drag coefficient 
		      betaP = 3./4.*rho*CD*magUr/ds*(deltaij+Kii*func_h*func_f);  	
                   
		    // calc particle's drag
                      drag = Volp*(betaP & Ur);

		    // Parcel approach: drag on one particle multiply by number of particles in a parcel
		    if(parcelApproach_)
		    {
			    drag *= npart_;
		    }	      
		      
                    if (modelType_=="B")
                        drag /= voidfraction;
                }

                if(verbose_ && index >=0 && index <2)
                {
                    Pout << "index = " << index << endl;
                    Pout << "Up = " << Up << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "rhos = " << rhos << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "Volp = " << Volp << endl;
                    Pout << "alps = " << alps << endl;
                    Pout << "CD = " << CD << endl;
		    Pout << "taup = " << taup << endl;
                    Pout << "delta_bar = " << delta_bar << endl;   
		    Pout << "func_f = " << func_f << endl;
                    Pout << "func_h = " << func_h << endl;
                    Pout << "Kii = " << Kii << endl;		    					
		    Pout << "betaP = " << betaP << endl;		   		    
                    Pout << "drag = " << drag << endl;
                }
            }
	    	    	    
            // set force on particle
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
	    	
              // 	if (modelType_=="B")
              // 		Cds()[index][0] = Volp*betaP/voidfraction;
              // 	else
              // 		Cds()[index][0] = Volp*betaP;	    	
            }
	    else
	    {
                for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];
            }
	    
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
