

#include <math.h>
#include "fvCFD.H"
#include "interpolationCellPoint.H"
#include "cfdemCloud.H"
#include "apiDensityModel.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//defineTypeNameAndDebug(apiDensityModel, 0);

//defineRunTimeSelectionTable(apiDensityModel, dictionary);

// * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

apiDensityModel::apiDensityModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm)
{
    //particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,maxCellsPerParticle_);
}


void apiDensityModel::apiDensityInterp( double** source, volScalarField& target )
{
   
   const int nparticles = this->particleCloud_.numberOfParticles();
   
   for( int i = 0; i < nparticles; ++i )
   {
	label cellI = this->particleCloud_.cellIDs()[i][0];
	if( cellI < 0 ) continue;
	
	scalar cellVolume = particleCloud_.mesh().V()[cellI];
	
	target[cellI] += source[i][0]/cellVolume;
	   
   } 
   
} 

}//end of namespace
