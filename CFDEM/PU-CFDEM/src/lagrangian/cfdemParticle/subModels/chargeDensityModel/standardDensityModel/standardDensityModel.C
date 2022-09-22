
#include "error.H"
#include "chargeDensityModel.H"
#include "standardDensityModel.H"
#include "locateModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(standardDensityModel, 0);

addToRunTimeSelectionTable
(
    chargeDensityModel,
    standardDensityModel,
    dictionary
);

standardDensityModel::standardDensityModel
(
    const dictionary& dict,
    cfdemCloud& sm
) : 
    chargeDensityModel( dict, sm )
{}


void standardDensityModel::setChargeDensity(double** const & pcharge_)
{    
    Info<<"Setting charge density... "<<particleCloud_.numberOfParticles()<<endl;
        
    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
    
        label cellI = particleCloud_.cellIDs()[index][0];
	
	if( cellI >= 0 )
	{
            scalar cellVolume = charge_density.mesh().V()[cellI];
	    scalar q = pcharge_[index][0];    
	    charge_density[cellI] += q/cellVolume;	    
	}
    
    }
    
}

































}
