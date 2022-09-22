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

#include "solidMassFluxController.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidMassFluxController, 0);

addToRunTimeSelectionTable
(
    forceModel,
    solidMassFluxController,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidMassFluxController::solidMassFluxController
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),      
    rhoParticle_(readScalar(propsDict_.lookup("rhoParticle"))),
    taup_(readScalar(propsDict_.lookup("taup"))),
    kappa_(readScalar(propsDict_.lookup("kappa"))),
    taud_(readScalar(propsDict_.lookup("taud"))),
    systemWeight_(readScalar(propsDict_.lookup("weight"))),
    initValue_( 0.0, 0.0, 0.0 ),
    momPartsOld_( 0,0,0 ),
    corDirection_( 0, 0, 1 ),
    coffA_( 0.1 ),
    corDirectionFlag_( false ),
    sourceReal_
    (   IOobject
        (
            "sourceReal",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0))  // N/m3 
    ),
    sourcePrev_
    (   IOobject
        (
            "sourcePrev",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0))  // N/m3 
    ),
    source_
    (   
        IOobject
        (
            "source",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	particleCloud_.mesh(),
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0)) 
    )      
{    
    if (propsDict_.found("verbose")) verbose_=true;        
    
    if (propsDict_.found("initialValue"))
    {
	initValue_ = vector( propsDict_.lookup("initialValue") );
	sourcePrev_.internalField() = initValue_;	
    }

    if (propsDict_.found("correctionDirection"))
    {
	corDirection_ = vector( propsDict_.lookup("correctionDirection") );
	corDirectionFlag_ = true;
	Info<<"Correction direction set on!"<<endl;
    }    
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

solidMassFluxController::~solidMassFluxController()
{}


scalar solidMassFluxController::controlFunction( vector momParts, vector momPartsPrev ) const
{
	
    scalar dMomParts = ( momParts - momPartsPrev ) & corDirection_;
    scalar coff( 0 );
    
    if( 
    	dMomParts * ( momParts & corDirection_ ) > 0 
       ) //change of momentum is in the same  direction as the momParts (source is increasing the momentum)
    {
        coff = 1.0;
    }else{ //change of momentum is in the opposite  direction as the momParts (source is decreasing the momentum)
        coff = 1.0-exp( -pow(dMomParts,2)/taud_ );
    }

    return coff;
	
}

tensor solidMassFluxController::controlFunctionT( vector momParts, vector momPartsPrev ) const
{
	
    tensor coff = tensor::zero;	
    scalar* coffE;
    
    		
    for( int i = 0; i < 3; ++i )
    {
    	
	vector e = vector::zero;
	
	if( i == 0 )
	{
	    
	    e.x() = 1.0;
	    coffE = &coff.xx();
	    
	}else if( i == 1 )
	{
	    
	    e.y() = 1.0;
	    coffE = &coff.yy();
	
	}else{
	    
	    e.z() = 1.0;
	    coffE = &coff.zz();
	
	}
	
	
	
	scalar dMomParts = ( momParts - momPartsPrev ) & e;

	if( 
    	    dMomParts * ( momParts & corDirection_ ) > 0 
	   ) //change of momentum is in the same  direction as the momParts (source is increasing the momentum)
	{
            *coffE = 1.0;
	}else{ //change of momentum is in the opposite  direction as the momParts (source is decreasing the momentum)
            *coffE = 1.0-exp( -pow(dMomParts,2)/taud_ );
	}
    
    }
    
    return coff;
	
} 

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidMassFluxController::setForce() const
{
    scalar ds(0);
    scalar mp(0);
    scalar Volp(0);
    scalar totVolp(0);    
    vector Up(0,0,0);
    vector momParts(0,0,0);
 
    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        ds   = 2*particleCloud_.radius(index);
	Volp = ds*ds*ds*M_PI/6;
	totVolp += Volp;
        mp = rhoParticle_*Volp;
	Up = particleCloud_.velocity(index);
	momParts += mp*Up; 
    }    	
    
    
    // Parallel computations
    reduce(totVolp , sumOp<scalar>());
    reduce(momParts  , sumOp<vector>());     
    
    if( corDirectionFlag_ )
    {
    	//source_.internalField() = ( ( sourcePrev_.internalField() + kappa_/taup_*momParts/totVolp ) & corDirection_ )
	//			* corDirection_ - ( taud_ * ( momParts/totVolp - momPartsOld_/totVolp ) & corDirection_ ) * corDirection_;
				
	/*source_.internalField() = ( sourcePrev_.internalField() & corDirection_ ) * corDirection_ +
				  (
				  ( kappa_/taup_*momParts/totVolp - taud_ * ( momParts/totVolp - momPartsOld_/totVolp )	) 
				  & corDirection_
				  ) * corDirection_;*/
				  
	vector a1 = momParts/totVolp;
	vector a2 = momPartsOld_/totVolp;			  
				  
	scalar coff = controlFunction( a1, a2 );			  
	source_.internalField() = ( sourcePrev_.internalField() & corDirection_ ) * corDirection_ +
				  coff * 	
				  (
				  	(  kappa_/taup_*momParts/totVolp ) & corDirection_
				  ) * corDirection_;
	
				  
	 //source_.internalField() = kappa_/taup_*momParts/totVolp;			   			  	
				
    }else{    
 
	vector a1 = momParts/totVolp;
	vector a2 = momPartsOld_/totVolp; 
        
	tensor coff = controlFunctionT( a1, a2 );
	
        source_.internalField() = sourcePrev_.internalField() +
				  ( coff & momParts * ( kappa_/taup_/totVolp ) );	

    }		
            	
    /*forAll(source_,cellI)
    {
        for(int j=0;j<3;j++ )
	{
	   if( source_[cellI][j] >= 0.5*systemWeight_ ) source_[cellI][j] = 0.5*systemWeight_;
	   if( source_[cellI][j] < -0.5*systemWeight_ ) source_[cellI][j] = -0.5*systemWeight_;
	}
    }*/
    	
	
    //source_ = -source_;				
    source_.correctBoundaryConditions();
    
    sourceReal_.internalField() = source_.internalField() + (-systemWeight_)*corDirection_;
    
    //JK added 9th July 2017
    particleCloud_.momCoupleM(1).setSourceField(sourceReal_);
    
    label cellI=0; ds = 0; Volp = 0 ;
        
    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle found on this processor
            {
                //Calc the particle volume
                ds   = 2*particleCloud_.radius(index);
                Volp = ds*ds*ds*M_PI/6;
		
                
		// set force on particle
		/*if(treatExplicit_) 
		{
	    	    for(int j=0;j<3;j++) expForces()[index][j] += Volp*source_[cellI][j];
        	}
		else
		{
	    	    for(int j=0;j<3;j++) impForces()[index][j] += Volp*source_[cellI][j];
		}*/

                // set force on particle
                for(int j=0;j<3;j++) 
                {
                    // calc particle's static pressure gradient force
		    DEMForces()[index][j] -= Volp*sourceReal_[cellI][j];
		    //DEMForces()[index][j] += Volp*source_[cellI][j];	
                }
		
            	if(verbose_ && index == 2)
		{
		  Info << "index = " << index << " Solid flux controller force ( "
                                                             << sourceReal_[cellI][0] << " "
                                                             << sourceReal_[cellI][1] << " "
                                                             << sourceReal_[cellI][2] << " )" 
							     << endl;
		}
		
            }else{
	    
	    	Info << "Warning: particle outside the processor!" <<endl;
		
		for(int j=0;j<3;j++)
		   DEMForces()[index][j] -= Volp*sourceReal_[0][j];
		   		
	    }
        //}
    }
    
    
    momPartsOld_.x() = momParts.x();
    momPartsOld_.y() = momParts.y();
    momPartsOld_.z() = momParts.z();
    
    // Update sourcePrev here 
    sourcePrev_ = source_;		
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
