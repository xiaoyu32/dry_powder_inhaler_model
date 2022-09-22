#ifndef EFIELD_MODEL_SCREENED_H
#define EFIELD_MODEL_SCREENED_H

#include "global_properties.h"	
#include "property_registry.h"
#include "math.h"
#include "stdlib.h"
	
#include "efield_model.h"

// -- class for electric field model --
class EfieldModelScreened : public EfieldModel
{
    
    protected :
    
    double screeningDistance;
        
    inline double gauss( const double r ) const
    {
        return exp( -r*r/ (2.0 * screeningDistance * screeningDistance) );
    }

    
    public :
    
    EfieldModelScreened() : 
    EfieldModel(),
    screeningDistance( 0 )
    {}

    ~EfieldModelScreened()
    {}


    virtual void connectToProperties(PropertyRegistry & registry) 
    {
       EfieldModel::connectToProperties( registry );
	
       // -- screening distance for P3M --
       registry.registerProperty("screeningDistance", &MODEL_PARAMS::createScreeningDistance);      
       registry.connect("screeningDistance", screeningDistance,"cohesion_model_electrostatic_screened");
	    
    }    

    virtual void computeEfieldCharge( double* normal, double r, double qi, double* E ) const
    {	
        //std::cout<<"EfieldModelScreened    permittivity is "<<permittivity<<std::endl;
    
	const double r2 = r*r;
	const double sqrt2sigma = sqrt( 2.0 ) * screeningDistance;
	const double sigma2 = 2.0 * screeningDistance * screeningDistance;
	
	const double erfcExp = erfc( r/sqrt2sigma );
	const double gaussExp = exp( -r2/sigma2 );

        const double Fn_coh = 1./(4.*M_PI*permittivity*r) *
		 		 ( 
				    erfcExp / r 
				  + sqrt( 2.0/3.1415926 ) * gaussExp / screeningDistance 
				 );    	
	
	//std::cout<<"Fn_coh = "<<Fn_coh<<" "<<"screeningDistance = "<<screeningDistance<<std::endl;
	
	for( int i = 0; i < 3; ++i )        
	{
	    E[i] += Fn_coh * qi * normal[i];
	}
	
    }
    
    virtual void computeGradEfieldCharge( double* normal, double r, double qi, double (&gradE)[3][3] ) const
    {
	const double r2 = r*r;
	const double sqrt2sigma = sqrt( 2.0 * screeningDistance );
	const double sigma2 = 2.0 * screeningDistance * screeningDistance;
	
	const double erfcExp = erfc( r/sqrt2sigma );
	const double gaussExp = exp( -r2/sigma2 );        
		
        const double Fn_coh = 1./(4.*M_PI*permittivity*r) *
		 		( 
				    erfcExp / r 
				  + sqrt( 2.0/M_PI ) * gaussExp / screeningDistance 
				 ); 
	
	const double gradFn_coh = 1./( 4.*M_PI*permittivity ) * 
				 ( 
				     sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * r * r ) 
				   + sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * screeningDistance * screeningDistance )
				   + erfcExp / ( r2 * r )
				 ); 

	// -- electric dipole of particle i --
 	const Tensor n( normal );
	
	const Tensor ntrans = n.transpose();
	
	const Tensor n2 = n * ntrans;
	
	Tensor efieldGrad = ( Tensor::eye() - n2 ) * qi * Fn_coh - n2 * qi * gradFn_coh;	
	
	// -- cpy tensor to the output array --
	for( int i = 0; i < 3; ++i )
	   for( int j = 0; j < 3; ++j )
	      gradE[i][j] += efieldGrad.get(i,j);
	
	
    }
    
    void computeEfieldPolarization( double* normal, double r, double* pi, double* E ) const
    {
        
	const double r2 = r*r;
	const double sqrt2sigma = sqrt( 2.0 * screeningDistance );
	const double sigma2 = 2.0 * screeningDistance * screeningDistance;
	
	const double erfcExp = erfc( r/sqrt2sigma );
	const double gaussExp = exp( -r2/sigma2 );
	
		
        const double Fn_coh = 1./(4.*M_PI*permittivity*r) *
		 		( 
				    erfcExp / r 
				  + sqrt( 2.0/M_PI ) * gaussExp / screeningDistance 
				 ); 
	
	const double gradFn_coh = 1./( 4.*M_PI*permittivity ) * 
				 ( 
				     sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * r * r ) 
				   + sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * screeningDistance * screeningDistance )
				   + erfcExp / ( r2 * r )
				 ); 

	// -- electric dipole of particle i --
	const Tensor pii( pi );
	
	//std::cout<<"pii = ";
	//pii.print();
	
 	const Tensor n( normal );
	
	//std::cout<<"n = ";
	//n.print();
	
	const Tensor ntrans = n.transpose();
	
	const Tensor n2 = n * ntrans;

	//std::cout<<"n2 = ";
	//n2.print();
	
	//std::cout<<"Fn_coh = "<<Fn_coh<<std::endl;
	//std::cout<<"gradFn_coh = "<<gradFn_coh<<std::endl;
	
	Tensor efield = ( Tensor::eye() - n2 ) * pii * Fn_coh + n2 * pii * gradFn_coh;
	
	//std::cout<<"efield = ";
	//efield.print();
	
	
	for( int i = 0; i < 3; ++i )
	{
	    E[i] += efield[i];
	}
	
    }
    
    void computeGradEfieldPolarization( double* normal, double r, double* pi, double (&gradE)[3][3] ) const
    {
        
	const double r2 = r*r;
	const double sqrt2sigma = sqrt( 2.0 * screeningDistance );
	const double sigma2 = 2.0 * screeningDistance * screeningDistance;
	
	const double erfcExp = erfc( r/sqrt2sigma );
	const double gaussExp = exp( -r2/sigma2 );
	
		
	const double gradFp_coh1 = 1.0/(4.*M_PI*permittivity*r) * 
				   (
				       erfcExp/ (r*r*r)
				     + sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * screeningDistance * screeningDistance)
				     + sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * r *r )
				   );
	
	
	const double gradFp_coh2 = 1.0/(4.*M_PI*permittivity) * 
		 		   (
				       3.0 * erfcExp/ (r2*r2)
				     + 3.0 * sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * r * r * r )
				     + sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * screeningDistance * screeningDistance * r )
				   );
	
	const double gradFp_coh3 = 1.0/(4.*M_PI*permittivity) * 
				   (
				       3.0 * erfcExp/(r2*r2) 
				     + 3.0 * sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * r * r * r)
				     + sqrt( 2.0/M_PI ) * gaussExp / ( screeningDistance * screeningDistance * screeningDistance * r )
				     + 3.0 * sqrt( 2.0/M_PI ) * gaussExp / ( pow(screeningDistance,4) )
				     + sqrt( 2.0/M_PI ) * gaussExp * r / ( pow(screeningDistance,5) )
				   );

	// -- electric dipole of particle i --
	const Tensor pii( pi );
 	const Tensor n( normal );
	
	const double pin = ( pii.transpose() * n ).value();
	
	const Tensor ntrans = n.transpose();
	
	const Tensor n2 = n * ntrans;

	const Tensor Psd1 = n2 * pin;

	const Tensor Psd2 = Tensor::eye() * pin 
		      	    + pii * ( n.transpose() ) 
		            - Psd1 * 2.0;

	const Tensor gradEfield =   Psd2 * gradFp_coh1
	 	      		   - ( pii * ( n.transpose() ) ) * gradFp_coh2
		      		   - Psd1 * gradFp_coh3 ;     
				   
	// -- cpy tensor to the output array --
	for( int i = 0; i < 3; ++i )
	   for( int j = 0; j < 3; ++j )
	      gradE[i][j] += gradEfield.get(i,j);				   
				       
    }
   
    
        

    // implicit coefficient for charge transfer model ( = 0 corresponds to explicit )
    double implCoeff( double radi, double radj )
    {
	const double sqrt2sigma = sqrt( 2.0 ) * screeningDistance;
	
	const double erfcExpi = erfc( radi/sqrt2sigma );
        const double erfcExpj = erfc( radj/sqrt2sigma );
	
        const double Fn_cohi = 1./(4.*M_PI*radi) *
		 		 ( 
				    erfcExpi / radi 
				  + sqrt( 2.0/M_PI ) * gauss(radi) / screeningDistance 
				 );    

        const double Fn_cohj = 1./(4.*M_PI*radj) *
		 		 ( 
				    erfcExpj / radi 
				  + sqrt( 2.0/M_PI ) * gauss(radj) / screeningDistance 
				 );    
     
       return Fn_cohi + Fn_cohj; 
       
    }
    
    // colliding particle electric field contribution for charge transfer model ( = 0 corresponds to explicit )
    double implElectricField( double qi, double qj, double radi, double radj )
    {

	const double sqrt2sigma = sqrt( 2.0 ) * screeningDistance;
	
	const double erfcExpi = erfc( radi/sqrt2sigma );
        const double erfcExpj = erfc( radj/sqrt2sigma );
	
        const double Fn_cohi = 1./(4.*3.1415926*permittivity*radi) *
		 		 ( 
				    erfcExpi / radi 
				  + sqrt( 2.0/M_PI ) * gauss(radi) / screeningDistance 
				 );    

        const double Fn_cohj = 1./(4.*3.1415926*permittivity*radj) *
		 		 ( 
				    erfcExpj / radi 
				  + sqrt( 2.0/M_PI ) * gauss(radj) / screeningDistance 
				 );
     
       return -qi * Fn_cohi + qj * Fn_cohj; 
    }
    

};


#endif
