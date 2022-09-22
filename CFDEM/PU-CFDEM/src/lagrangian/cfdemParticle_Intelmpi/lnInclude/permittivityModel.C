

#include <math.h>
#include "fvCFD.H"
#include "interpolationCellPoint.H"
#include "cfdemCloud.H"
#include "permittivityModel.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PermittivityModel, 0);

defineRunTimeSelectionTable(PermittivityModel, dictionary);

// * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




PermittivityModel::PermittivityModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    permittivitySolid_
    (
	dict.subDict
	(
	  "permittivityProps"
	).lookupOrDefault("solid",1.)    
    ),
    permittivityGas_
    (
	dict.subDict
	(
	  "permittivityProps"
	).lookupOrDefault("gas",1.)    
    )
{
    //particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,maxCellsPerParticle_);
}

void PermittivityModel::setPermittivity( const volScalarField& alphap, volScalarField& permittivity ) const
{
    
    forAll( alphap, iCell )
    {
        permittivity[iCell] = calcPermittivity( alphap[iCell] );
    }
    
}

double PermittivityModel::calcPermittivity( const double alphap ) const
{
    return ( alphap * permittivitySolid_ + (1-alphap) * permittivityGas_ );
}


}//end of namespace
