
#ifndef chargeDensityModel_H
#define chargeDensityModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

namespace Foam
{


class chargeDensityModel
{
    
protected:
    
    const dictionary& dict_ ;
    cfdemCloud& particleCloud_;
         
    mutable volScalarField charge_density; 

public:
    
    TypeName("chargeDensityModel");
    
    declareRunTimeSelectionTable
    (
       autoPtr,
       chargeDensityModel,
       dictionary,
       (
           const dictionary& dict,
           cfdemCloud& sm
        ),
       (dict,sm)
    );
   
   
       // Constructors

        //- Construct from components
        chargeDensityModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

    virtual ~chargeDensityModel(){};


    // Selector

    static autoPtr<chargeDensityModel> New
    (
       const dictionary& dict,
       cfdemCloud& sm
    );
    

    inline double n3root( double );  
    void compEfCorrection( double** & ef_cpr ); // deprecated method (used in earlier version of P3M)
    void getElectricField( double** & pef_, volVectorField& ef  ); // --- should be moved somewhere else --- (has nothing to do with charge density)

    virtual void setChargeDensity(double ** const &) = 0;
    void resetChargeDensity();
    tmp<volScalarField> chargeDensityInterp() const;
   
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
