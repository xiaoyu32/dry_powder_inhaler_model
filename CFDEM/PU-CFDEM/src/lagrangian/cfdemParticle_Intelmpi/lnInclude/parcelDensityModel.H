#ifndef parcelDensityModel_H
#define parcelDensityModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "chargeDensityModel.H"

namespace Foam
{


class parcelDensityModel : public chargeDensityModel
{
    
protected:
    
    int apiType_;
    int nParcel_;

public:
    
    TypeName("parcel");

       // Constructors

        //- Construct from components
        parcelDensityModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

    virtual ~parcelDensityModel(){};

    virtual void setChargeDensity(double ** const &);
   
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
















#endif
