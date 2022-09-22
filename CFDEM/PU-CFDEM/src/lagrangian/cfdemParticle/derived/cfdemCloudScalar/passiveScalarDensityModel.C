

#include <math.h>
#include "fvCFD.H"
#include "interpolationCellPoint.H"
#include "cfdemCloud.H"
#include "passiveScalarDensityModel.H"
#include "List.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(passiveScalarDensityModel, 0);

defineRunTimeSelectionTable(passiveScalarDensityModel, dictionary);

// * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




passiveScalarDensityModel::passiveScalarDensityModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    kDissolve_(NULL),
    maxCellsPerParticle_(1),
    kHenry
    ( 
        dict.subDict("passiveScalarTransferModel").lookup("henrysConstant")
    ),
    Tcarrier
    (
        readScalar(dict.subDict("passiveScalarTransferModel").lookup("carrierGasTemperature"))
    ),
    SchmidtNumber
    (
        readScalar(dict.subDict("passiveScalarTransferModel").lookup("SchmidtNumber"))
    ),
    viscosity
    (
        readScalar(dict.subDict("passiveScalarTransferModel").lookup("viscosity"))
    )      
{
    particleCloud_.dataExchangeM().allocateArray(kDissolve_,1,maxCellsPerParticle_);
}


passiveScalarDensityModel::~passiveScalarDensityModel()
{
    particleCloud_.dataExchangeM().destroy(kDissolve_,1);
}



inline double Foam::passiveScalarDensityModel::n3root( double a )
{
    return std::pow(a, 1.0/3.0 );
}


void Foam::passiveScalarDensityModel::getScalarTransfer
(  
	List<double**> & pscalar_, 
	List< autoPtr<volScalarField> >& expPassiveScalarSource, 
	List< autoPtr<volScalarField> >& impPassiveScalarSource,
	volScalarField& alpha,
	volVectorField& Us,
	volVectorField& U  
)
{
    
    reAllocArrays();
    
    //reset particle volume in cell
    //particleVolumeInCell.internalField() = 0.0;
    
    for( label ii=0; ii < numberOfPassiveScalars(); ++ii )
    { 
	expPassiveScalarSource[ii]().internalField() = 0.0;
	impPassiveScalarSource[ii]().internalField() = 0.0;
    }
    
    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
    
        label cellI = particleCloud_.cellIDs()[index][0];
	
	if( cellI > 0 )
	{
	     
	     double pradii = particleCloud_.radii()[index][0];
	     
	     //cell volume
	     const double Vcell = expPassiveScalarSource[0]().mesh().V()[cellI];
	     
	     //particle area
	     const double Api = 4.0 * 3.1415926 * pradii * pradii;
	     
	     //time step
	     const double dt = expPassiveScalarSource[0]().mesh().time().deltaT().value();
	     
	     //compute interfacial transfer coefficient
	     //const double kDissolve = updateInterfacialTransferCoefficients( alpha[cellI], Us[cellI], U[cellI], 2.0*pradii );
	     const vector vp = particleCloud_.velocity(index);
	     const double kDissolve = updateInterfacialTransferCoefficients( alpha[cellI], vp, U[cellI], 2.0*pradii );
	     
	     //save the interfacial transfer coefficient for computing the flux to particles
	     kDissolve_[index][0] = kDissolve;
	     
	     //loop over passive scalars
	     for( label ii = 0; ii < numberOfPassiveScalars(); ++ii )
	     {
	         
		 const double chi = 1.0/(Rgas * Tcarrier * kHenry[ii]);
		 
		 const double imp_coff = 1.0 + 3.0/pradii * kDissolve* chi * dt;

		 //explicit flux from partilces
		 const double cp_exp = 1/Vcell * Api * kDissolve * chi * pscalar_[ii][index][0]/imp_coff;

		 //implicit flux (coefficient of the fluid passive scalar concentration
		 //FIXED IT: changed the sign in the implicit coefficient 1 + ... to 1 - ... (JK)
		 const double cg_imp = -1/Vcell * Api * kDissolve *  
	     			       (1 - 3/pradii * kDissolve * chi * dt / imp_coff);
				       
		 expPassiveScalarSource[ii]()[cellI] += cp_exp;
		 impPassiveScalarSource[ii]()[cellI] += cg_imp;  
		 
	     }  
	}else
	{
	     kDissolve_[index][0] = 0.0;
	}
	    
    }
    
    
}


void Foam::passiveScalarDensityModel::computeScalarTransfer
( 
	List<double**> & pscalar_, 
	List<double**> & dpscalar_, 
	List< autoPtr<volScalarField> >& passiveScalar 
)
{
    
    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
        for( label ii = 0; ii < numberOfPassiveScalars(); ++ii ) dpscalar_[ii][index][0] = 0.0;
    }	    
    
    
    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
        
	label cellI = particleCloud_.cellIDs()[index][0];
	
	double pradii = particleCloud_.radii()[index][0];
	
	const double dt = passiveScalar[0]().mesh().time().deltaT().value();
	
	if( cellI >= 0 )
	{
	    
	    const double kDissolve =  kDissolve_[index][0];
	    
	    for( label ii = 0; ii < numberOfPassiveScalars(); ++ii )
	    {
	    	
		const double chi = 1.0/(Rgas * Tcarrier * kHenry[ii]);
		
		const double imp_coff = 1.0 + 3.0/pradii * kDissolve * chi * dt;
		const double cg = passiveScalar[ii]()[cellI];
		const double cp_prev = pscalar_[ii][index][0];

		const double cp_next = ( cp_prev + 3/pradii * kDissolve * dt * cg )/imp_coff;

		dpscalar_[ii][index][0] += 3.0/pradii * kDissolve * ( 
	    				     		cg - chi * cp_next 
								    );
					
	    }
	}
	
    }
    
}

double Foam::passiveScalarDensityModel::updateInterfacialTransferCoefficients 
( 
     	const double& vof,
	const vector& Us,
	const vector& U,
	const double dp
) const
{
    
    //local Reynolds number
    const scalar Re = vof * mag( U - Us ) * dp / viscosity;
    
    //particles Sherwood number
    const scalar SherwoodNumber = (7.0 - 10.0 * vof + 5.0 * vof * vof ) * 
    				  ( 1.0 + 0.7 * pow(Re,0.2) * pow( SchmidtNumber, 1.0/3.0) ) + 
				  (1.33 - 2.4 * vof + 1.3 * vof * vof ) *  
				  pow(Re,0.7) * pow( SchmidtNumber, 1.0/3.0);
    
    const double massDiffusivity = viscosity/SchmidtNumber;
     
    //(convective) mass transfer coefficient
    const double coeff = SherwoodNumber * massDiffusivity / dp;
    
    return coeff;    
    
}


int Foam::passiveScalarDensityModel::numberOfPassiveScalars()
{
    return (int)kHenry.size();
}


void Foam::passiveScalarDensityModel::reAllocArrays() const
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        // get arrays of new length
        particleCloud_.dataExchangeM().allocateArray(kDissolve_,0,1);
    }
}

}//end of namespace
