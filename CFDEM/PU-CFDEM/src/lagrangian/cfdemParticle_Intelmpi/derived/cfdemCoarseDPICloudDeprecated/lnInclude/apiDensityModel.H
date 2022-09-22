
#ifndef apiDensityModel_H
#define apiDensityModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

namespace Foam
{


class apiDensityModel
{
    
protected:
    
    const dictionary& dict_ ;
    cfdemCloud& particleCloud_;
        
    mutable double ** cellsPerParticle_;
    int maxCellsPerParticle_;

public:
    
    /*TypeName("apiDensityModel");
    
    declareRunTimeSelectionTable
    (
       autoPtr,
       apiDensityModel,
       dictionary,
       (
           const dictionary& dict,
           cfdemCloud& sm
        ),
       (dict,sm)
    );*/
   
   
       // Constructors

        //- Construct from components
        apiDensityModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

    virtual ~apiDensityModel(){};


    // Selector

    static autoPtr<apiDensityModel> New
    (
       const dictionary& dict,
       cfdemCloud& sm
    );
    
    void apiDensityInterp( double**, volScalarField& );
   
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
