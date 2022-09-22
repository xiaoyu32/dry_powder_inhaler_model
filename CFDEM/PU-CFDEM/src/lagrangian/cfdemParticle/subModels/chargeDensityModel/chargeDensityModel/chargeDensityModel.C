

#include <math.h>
#include "fvCFD.H"
#include "interpolationCellPoint.H"
#include "cfdemCloud.H"
#include "chargeDensityModel.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(chargeDensityModel, 0);

defineRunTimeSelectionTable(chargeDensityModel, dictionary);

// * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




chargeDensityModel::chargeDensityModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    charge_density
    (   IOobject
        (
            "rhoe",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volScalarField> ("rhoe")
        /*sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 1)*/
    )
{
    //particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,maxCellsPerParticle_);
}

inline double chargeDensityModel::n3root( double a )
{
    return std::pow(a, 1.0/3.0 );
}


void chargeDensityModel::compEfCorrection( double** & ef_cpr )
{
    
         
    int i,j,k,n;
        
    //allocate
    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
	
	for( k = 0; k < 3; ++k )
	{
	    ef_cpr[index][k] = 0.0;
	}

        //for( k = 0; k < 37; ++k )
	//{
	//    cell_l[index][k] = 0.0;
	//}
    } 
    
    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {	    
        
	label cellI = particleCloud_.cellIDs()[index][0];
	
	//particle is on the processor
	if( cellI >= 0 )
	{
	    //cell_l[index][0] = n3root( charge_density.mesh().V()[cellI] ) * n3root( 0.238732343 );	
            
	    /*const cell& faces = charge_density.mesh().cells()[cellI];	
	    
            //hexahedra surface areas and face center points
	    k = 0;
	    for( j = 0; j < 6; ++j )
	    {

  		for( int kk = 0; kk < 3; ++kk )
		{ 
		    cell_l[index][k] = charge_density.mesh().Cf()[faces[j]][kk];
		    ++k;
            	}
		
		for( int kk = 0; kk < 3; ++kk )
		{
		    cell_l[index][k] = charge_density.mesh().Sf()[faces[j]][kk];
		    ++k;
		}

	    }*/
	    
	    
            //cell volume
            //cell_l[index][k] = charge_density.mesh().V()[cellI];	    

	    	
	    //save particles cell centers
	    for( j = 0; j < 3; ++j )
	    {
	        ef_cpr[index][j] = charge_density.mesh().C()[cellI][j];
	    } 	     
	    
	}

    }
     
}

void chargeDensityModel::getElectricField( double** & pef_, volVectorField& electric_field  )
{
    
    vector position(0,0,0);
    vector ef(0,0,0);
    //interpolationCellPoint<vector> efInterpolator_( electric_field );

    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
        for( int j = 0; j < 3; ++j )
        {
 	    pef_[index][j] = 0.0;
        }
    }		

    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {	    
        label cellI = particleCloud_.cellIDs()[index][0];
	ef = vector(0,0,0);

	if( cellI >= 0 )
	{
            position = particleCloud_.position(index); 	     
	    //ef = efInterpolator_.interpolate(position,cellI); //interpolation
	    
	    for( int j = 0; j < 3; ++j )
	    {
                ef[j] = electric_field[cellI][j]; //injection
            }	

            for( int j = 0; j < 3; ++j )
	    {
 		pef_[index][j] = ef[j];
      	    }
 
	}

    }


}

/* JK: moved to child class
void ChargeDensityModel::setChargeDensity(double** const & pcharge_)
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
*/

void Foam::chargeDensityModel::resetChargeDensity()
{
    //set zero charge density 
    
    charge_density.internalField() = 0.0;
}

tmp<volScalarField> Foam::chargeDensityModel::chargeDensityInterp() const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "rhoe",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimensionSet(0, 0, 0, 0, 0),
                0
            )
        )
    );
    
    tsource() = charge_density;

    return tsource;

}

}//end of namespace
