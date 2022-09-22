#ifndef EFIELD_MODEL_NORMAL_H
#define EFIELD_MODEL_NORMAL_H

#include "global_properties.h"	
#include "efield_model.h"

// -- class for electric field model --
class EfieldModelNormal : public EfieldModel
{
    
    protected :
    
    public :
    
    EfieldModelNormal() : 
    EfieldModel()
    {}
    
    ~EfieldModelNormal()
    {}

    virtual void computeEfieldCharge( double* normal, double r, double qi, double* E ) const
    {
        //std::cout<<"permittivity is "<<permittivity<<std::endl;
	
	const double cof = 1.0/(4*M_PI*permittivity);
    
	for( int i = 0; i < 3; ++i )        
	{
	    E[i] += cof * qi * normal[i] / (r*r);
	}
	
    }

    virtual void computeGradEfieldCharge( double* normal, double r, double qi, double (&gradE)[3][3] ) const
    {
        const double cof = 1.0/(4*M_PI*permittivity);
 	
	const Tensor n( normal );
	const Tensor ntrans = n.transpose();
	const Tensor n2 = n * ntrans;	
	
	Tensor efieldGrad = ( n2 * 3.0 -  Tensor::eye() ) * ( qi * cof / (r*r*r) );
	
	// -- cpy tensor to the output array --
	for( int i = 0; i < 3; ++i )
	   for( int j = 0; j < 3; ++j )
	      gradE[i][j] += efieldGrad.get(i,j);
	      	
    }

    void computeEfieldPolarization( double* normal, double r, double* pi, double* E ) const
    {
        
	const double cof = 1.0/(4*M_PI*permittivity);
	
	double pin = 0;
	 
	for( int i = 0; i < 3; ++i )
	{
	   pin += pi[i] * normal[i];
	} 
	
	for( int i = 0; i < 3; ++i )
	{
	    E[i] += cof * ( 3.0 * pin * normal[i] - pi[i] )/(r*r*r);
	}
	
    }

    void computeGradEfieldPolarization( double* n, double r, double* p, double (&gradE)[3][3] ) const
    {
        const double cof = 1.0/(4.0*M_PI*permittivity);   
	
	double pin = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
	
	for( int i = 0; i < 3; ++i )
	{
	    for( int j = 0; j < 3; ++j )
	    {
		gradE[i][j] += cof * (  
			       		 3*n[i]*p[j]/(r*r*r*r)
			     	       + 3*p[i]*n[j]/(r*r*r*r)
			               - 15*pin*n[i]*n[j]/(r*r*r*r) 
			             );	    
	       // -- isotropic part --		   
	       if( i == j )
	       {
		   gradE[i][j] += cof * 3 * pin/(r*r*r*r);
	       }			   
			   
	    }
	}    
	
    }

    // implicit coefficient for charge transfer model ( = 0 corresponds to explicit )
    double implCoeff( double radi, double radj )
    { 
       return 1.0/(4*M_PI*radi*radi) + 1.0/(4*M_PI*radj*radj); 
    }
    
    // colliding particle electric field contribution for charge transfer model ( = 0 corresponds to explicit )
    double implElectricField( double qi, double qj, double radi, double radj )
    { 
       return -qi/(4*M_PI*permittivity*radi*radi) + qj/(4*M_PI*permittivity*radj*radj); 
    }


};


#endif
