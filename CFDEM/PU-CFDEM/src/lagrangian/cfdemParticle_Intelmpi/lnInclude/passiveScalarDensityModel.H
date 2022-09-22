
#ifndef PassiveScalarDensityModel_H
#define PassiveScalarDensityModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "List.H"


namespace Foam
{


class passiveScalarDensityModel
{
    
protected:
    
    const dictionary& dict_ ;
    cfdemCloud& particleCloud_;
   
    mutable double ** kDissolve_;
    
    int maxCellsPerParticle_;
    
    //temperature
    mutable double Tcarrier;
    
    //universal gas constant
    static const double Rgas = 8.3144598;
    
    const double viscosity;
    const double SchmidtNumber;
    
    //Henry's constant for solubility
    //mutable double kHenry;
    mutable scalarList kHenry;

    double updateInterfacialTransferCoefficients
    ( 
     	    const double& alpha,
	    const vector& Us,
	    const vector& U,
	    const double dp
    ) const;

public:
    
    TypeName("passiveScalarDensityModel");
    
    declareRunTimeSelectionTable
    (
       autoPtr,
       passiveScalarDensityModel,
       dictionary,
       (
           const dictionary& dict,
           cfdemCloud& sm
        ),
       (dict,sm)
    );
   
   
       // Constructors

        //- Construct from components
        passiveScalarDensityModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

    virtual ~passiveScalarDensityModel();


    // Selector

    static autoPtr<passiveScalarDensityModel> New
    (
       const dictionary& dict,
       cfdemCloud& sm
    );
    

    inline double n3root( double );  
    
    //explicit and implicit parts
    void getScalarTransfer( 
    				List<double**> & pscalar_, 
				List< autoPtr<volScalarField> >& expPassiveScalarSource, 
				List< autoPtr<volScalarField> >& impPassiveScalarSource,
				volScalarField& alpha,
				volVectorField& Us,
				volVectorField& U 
			);
    
    int numberOfPassiveScalars();
    
    //
    void computeScalarTransfer( 
    				List<double**> &, 
				List<double**> &, 
				List< autoPtr<volScalarField> >& 
				);
    
    int maxCellsPerParticle() const;
    void reAllocArrays() const;    
   
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
