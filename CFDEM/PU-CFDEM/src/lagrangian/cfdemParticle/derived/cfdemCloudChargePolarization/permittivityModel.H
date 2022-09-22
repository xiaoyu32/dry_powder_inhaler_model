
#ifndef PermittivityModel_H
#define PermittivityModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

namespace Foam
{


class PermittivityModel
{
    
protected:
    
    const dictionary& dict_ ;
    const cfdemCloud& particleCloud_;

    const double permittivitySolid_;
    const double permittivityGas_;

public:
    
    TypeName("PermittivityModel");
    
    declareRunTimeSelectionTable
    (
       autoPtr,
       PermittivityModel,
       dictionary,
       (
           const dictionary& dict,
           cfdemCloud& sm
        ),
       (dict,sm)
    );
   
   
       // Constructors

        //- Construct from components
        PermittivityModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

    virtual ~PermittivityModel(){};


    // Selector

    static autoPtr<PermittivityModel> New
    (
       const dictionary& dict,
       cfdemCloud& sm
    );
    

    void setPermittivity( const volScalarField& alphap, volScalarField& permittivity ) const;
    
    protected :
    virtual double calcPermittivity( const double alphap ) const;
    
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
