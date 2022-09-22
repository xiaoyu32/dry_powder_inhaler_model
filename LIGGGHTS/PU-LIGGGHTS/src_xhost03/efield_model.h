#ifndef EFIELD_MODEL_H
#define EFIELD_MODEL_H

// -- macros --
#define ADDEFIELDMODEL(name,type,efieldModel) if( strcmp(arg[iarg],#name) == 0 ){ efieldModel=new type(); std::cout<<"Using model: "<<#name<<std::endl; }	

#include "global_properties.h"				
#include "property_registry.h"				
#include "math.h"

					
class EfieldModelNormal;
class EfieldModelScreened;

// #define M_PI (3.1415926)

// -- class for electric field model --
class EfieldModel
{
    
    protected :
    
    // -- class for tensor calculations --
    #include "tensorOperations.h"
    
    double permittivity;
    
    
    public :
    
    EfieldModel() :
    permittivity( 0 )
    {}

    ~EfieldModel()
    {}
    
    virtual void connectToProperties(PropertyRegistry & registry) 
    {
       registry.registerProperty("permittivity", &MODEL_PARAMS::createPermittivity);      
       registry.connect("permittivity", permittivity,"efieldModel");
       // registry.connect("screeningDistance", screeningDistance,"cohesion_model_electrostatic_screened");        
    }    
    
    void setPermittivity( double permittivity_ )
    {
        permittivity = permittivity_;
    }
    
    virtual void computeEfieldCharge( double* n, double r, double qi, double* E ) const
    {}

    virtual void computeGradEfieldCharge( double* n, double r, double qi, double (&gradE)[3][3] ) const
    {}

    virtual void computeEfieldPolarization( double* n, double r, double* pi, double* E ) const
    {}

    virtual void computeGradEfieldPolarization( double* n, double r, double* pi, double (&gradE)[3][3] ) const
    {}
    
    // implicit coefficient for charge transfer model ( = 0 corresponds to explicit )
    virtual double implCoeff( double radi, double radj )
    { return 1.0; }
    
    // colliding particle electric field contribution for charge transfer model ( = 0 corresponds to explicit )
    virtual double implElectricField( double qi, double qj, double radi, double radj )
    { return 0; }
    
        
};


#endif












