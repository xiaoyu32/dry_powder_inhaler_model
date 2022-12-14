
#ifndef chargeDensityModel_H
#define chargeDensityModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

namespace Foam
{


class ChargeDensityModel
{
    
protected:
    
    const dictionary& dict_ ;
    cfdemCloud& particleCloud_;
         
    mutable volScalarField charge_density; 

    mutable double ** cellsPerParticle_;
    int maxCellsPerParticle_;

public:
    
    TypeName("ChargeDensityModel");
    
    declareRunTimeSelectionTable
    (
       autoPtr,
       ChargeDensityModel,
       dictionary,
       (
           const dictionary& dict,
           cfdemCloud& sm
        ),
       (dict,sm)
    );
   
   
       // Constructors

        //- Construct from components
        ChargeDensityModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

    virtual ~ChargeDensityModel(){};


    // Selector

    static autoPtr<ChargeDensityModel> New
    (
       const dictionary& dict,
       cfdemCloud& sm
    );
    

    inline double n3root( double );  
    void compEfCorrection( double** & ef_cpr );
    void getElectricField( double** & pef_, volVectorField& ef  );

    void setChargeDensity(double ** const &);
    void resetChargeDensity();
    int maxCellsPerParticle() const;
    void reAllocArrays() const;    
    tmp<volScalarField> chargeDensityInterp() const;
   
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
