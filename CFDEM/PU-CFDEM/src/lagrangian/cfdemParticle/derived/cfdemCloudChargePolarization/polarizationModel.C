

#include <math.h>
#include "fvCFD.H"
#include "interpolationCellPoint.H"
#include "cfdemCloud.H"
#include "polarizationModel.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PolarizationModel, 0);

defineRunTimeSelectionTable(PolarizationModel, dictionary);

// * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




PolarizationModel::PolarizationModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    cellsPerParticle_(NULL),
    maxCellsPerParticle_(1)
{
    //particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,maxCellsPerParticle_);
}

inline double PolarizationModel::n3root( double a )
{
    return std::pow(a, 1.0/3.0 );
}


void PolarizationModel::getElectricFieldGradient( double** & efGrad_, volTensorField& gradElectricField )
{
    
    vector position(0,0,0);
    //interpolationCellPoint<vector> efInterpolator_( electric_field );

    reAllocArrays();

    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
        for( int j = 0; j < 9; ++j )
        {
 	    efGrad_[index][j] = 0.0;
        }
    }		

    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {	    
        label cellI = particleCloud_.cellIDs()[index][0];

	if( cellI >= 0 )
	{
            position = particleCloud_.position(index); 	     
	    //ef = efInterpolator_.interpolate(position,cellI); //interpolation
	    

            efGrad_[index][0] = gradElectricField[cellI].xx(); //injection
 	    efGrad_[index][1] = gradElectricField[cellI].xy();
	    efGrad_[index][2] = gradElectricField[cellI].xz();
	    
	    efGrad_[index][3] = gradElectricField[cellI].yx();
	    efGrad_[index][4] = gradElectricField[cellI].yy();
	    efGrad_[index][5] = gradElectricField[cellI].yz();

	    efGrad_[index][6] = gradElectricField[cellI].zx();
	    efGrad_[index][7] = gradElectricField[cellI].zy();
	    efGrad_[index][8] = gradElectricField[cellI].zz();
	    
	}

    }


}

void Foam::PolarizationModel::reAllocArrays() const
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        // get arrays of new length
        particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,1);
    }
}

int Foam::PolarizationModel::maxCellsPerParticle() const
{
    return maxCellsPerParticle_;
}

}//end of namespace
