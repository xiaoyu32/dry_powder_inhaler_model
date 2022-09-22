
#ifndef PolarizationModel_H
#define PolarizationModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"



namespace Foam
{


class PolarizationModel
{
    
protected:
    
    const dictionary& dict_ ;
    cfdemCloud& particleCloud_;
         
    mutable double ** cellsPerParticle_;
    int maxCellsPerParticle_;

public:
    
    TypeName("PolarizationModel");
    
    declareRunTimeSelectionTable
    (
       autoPtr,
       PolarizationModel,
       dictionary,
       (
           const dictionary& dict,
           cfdemCloud& sm
        ),
       (dict,sm)
    );
   
   
       // Constructors

        //- Construct from components
        PolarizationModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

    virtual ~PolarizationModel(){};


    // Selector

    static autoPtr<PolarizationModel> New
    (
       const dictionary& dict,
       cfdemCloud& sm
    );
    

    inline double n3root( double );  

    int maxCellsPerParticle() const;
    void reAllocArrays() const;    
    void getElectricFieldGradient( double** & efGrad_, volTensorField& gradElectricField );
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
