#ifndef standardDensityModel_H
#define standardDensityModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "chargeDensityModel.H"

namespace Foam
{


class standardDensityModel : public chargeDensityModel
{
    
public:
    
    TypeName("standard");

       // Constructors

        //- Construct from components
        standardDensityModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

    virtual ~standardDensityModel(){};

    virtual void setChargeDensity(double ** const &);
   
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
















#endif
