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
#include "apiTransferModel.H"
#include "vanWachem.H"
#include "cfdemCloud.H"
#include "cfdemCoarseDPICloud.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vanWachemTransferModel, 0);

addToRunTimeSelectionTable
(
    apiTransferModel,
    vanWachemTransferModel,
    dictionary
);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
vanWachemTransferModel::vanWachemTransferModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
apiTransferModel( dict, sm ),
E0
(
   readScalar(dict.lookup("E0"))
),
E1
(
   readScalar(dict.lookup("E1"))
),
U0
(
   readScalar(dict.lookup("U0"))
)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

vanWachemTransferModel::~vanWachemTransferModel()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

double vanWachemTransferModel::diffNorm( const double* a, const double* b) const
{
    
    double res = 0;
    for( int i = 0; i < 3; ++i )
       res += pow( a[i]-b[i], 2 );
    return sqrt( res );
}

// -- explicit API particle source -- (+ apiExplFlux_)
void vanWachemTransferModel::explicitSources( volScalarField& explicitFlux )
{

    const int nparticles = this->particleCloud_.numberOfParticles();
    double** particleVel = this->particleCloud_.velocities();
    double** fluidVel = this->particleCloud_.fluidVels();
    double** radii = this->particleCloud_.radii();
    
    cfdemCoarseDPICloud& myCloud = dynamic_cast<cfdemCoarseDPICloud&>( particleCloud_ );
    
    double** carrierC = myCloud.carrierConcentration();
    
    for( int i = 0; i < nparticles; ++i )
    {
        
	label cellI = this->particleCloud_.cellIDs()[i][0];
	if( cellI < 0 ) continue;
	
	// -- carrier particle surface area --
	scalar Api = 4.0 * M_PI * pow( radii[i][0], 2 );
	
	double vel
	(
	   diffNorm(  particleVel[i], fluidVel[i] )
	);
	
	scalar ff
	(
	   (vel-U0) > 0 ?  vel-U0 : 0
	);
	
	scalar cellVolume = particleCloud_.mesh().V()[cellI];
	
	if( ff > 0 )
	{
	    explicitFlux[cellI] += E0 * Foam::pow(ff,E1) * Api * carrierC[i][0] / cellVolume;
	}
	
    }
    
}

// -- coefficient of fluid api particle concentration in the transport equation ( - fvm::Sp( apiImplFlux_, apiFluidDensity_ ) --
void vanWachemTransferModel::implicitSources( volScalarField& implicitFlux )
{
    const int nparticles = this->particleCloud_.numberOfParticles();
    double** particleVel = this->particleCloud_.velocities();
    double** fluidVel = this->particleCloud_.fluidVels();
    double** radii = this->particleCloud_.radii();
    double** particleVoidfractions = this->particleCloud_.voidfractions();
        
    for( int i = 0; i < nparticles; ++i )
    {
        
	label cellI = this->particleCloud_.cellIDs()[i][0];
	if( cellI < 0 ) continue;
	
	// -- carrier particle surface area --
	scalar Api = 4.0 * M_PI * pow( radii[i][0], 2 );
	scalar cellVolume = particleCloud_.mesh().V()[cellI];
	
	scalar vel
	(
	   diffNorm(  particleVel[i], fluidVel[i] )
	);

	implicitFlux[cellI] += vel * Api * particleVoidfractions[i][0] / cellVolume;
	
    }
    
}

// -- compute number of API particles transferred to each carrier particle per unit time --
void vanWachemTransferModel::computeSource( const volScalarField& apiDensity, double** carrierFlux )
{
    const int nparticles = this->particleCloud_.numberOfParticles();
    double** particleVel = this->particleCloud_.velocities();
    double** fluidVel = this->particleCloud_.fluidVels();
    double** radii = this->particleCloud_.radii();
    double** particleVoidfractions = this->particleCloud_.voidfractions();
 
    cfdemCoarseDPICloud& myCloud = dynamic_cast<cfdemCoarseDPICloud&>( particleCloud_ );   
    double** carrierC = myCloud.carrierConcentration(); 
        
    for( int i = 0; i < nparticles; ++i )
    {
        
	// -- implicit source --
	label cellI = this->particleCloud_.cellIDs()[i][0];
	if( cellI < 0 ) continue;
	
	// -- carrier particle surface area --
	scalar Api = 4.0 * M_PI * pow( radii[i][0], 2 );
	
	scalar vel
	(
	   diffNorm(  particleVel[i], fluidVel[i] )
	);

	scalar implSource = vel * Api * particleVoidfractions[i][0] * apiDensity[cellI];	
	
	// -- explicit source --
	scalar ff
	(
	   (vel-U0) > 0 ?  vel-U0 : 0
	);
	
	scalar cellVolume = particleCloud_.mesh().V()[cellI];
	scalar explSource = -E0 * Foam::pow(ff,E1) * Api * carrierC[i][0];
	
	carrierFlux[i][0] = implSource + explSource;
		
    }
     
}

} // End namespace Foam

// ************************************************************************* //









