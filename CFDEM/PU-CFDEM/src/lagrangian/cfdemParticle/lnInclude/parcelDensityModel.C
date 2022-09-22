
#include "error.H"
#include "chargeDensityModel.H"
#include "parcelDensityModel.H"
#include "locateModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(parcelDensityModel, 0);

addToRunTimeSelectionTable
(
    chargeDensityModel,
    parcelDensityModel,
    dictionary
);

parcelDensityModel::parcelDensityModel
(
    const dictionary& dict,
    cfdemCloud& sm
) : 
    chargeDensityModel( dict, sm ),
    nParcel_
    ( 
  	dict.subDict(typeName + "Props").lookupOrDefault("nParcel", 1.0)
    ),
    apiType_
    (
        dict.subDict(typeName + "Props").lookupOrDefault("apiType", 1)
    )
{}


void parcelDensityModel::setChargeDensity(double** const & pcharge_)
{    
    Info<<"Setting charge density... "<<particleCloud_.numberOfParticles()<<endl;
        
    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
    
        label cellI = particleCloud_.cellIDs()[index][0];
	
	if( cellI >= 0 )
	{
	
	    int ptype = particleCloud_.types()[index][0];
	
            scalar cellVolume = charge_density.mesh().V()[cellI];
	    scalar q = ptype == apiType_  ? nParcel_ * pcharge_[index][0] : pcharge_[index][0];
	    
	    charge_density[cellI] += q/cellVolume;	    
	}
    
    }
    
}

































}
