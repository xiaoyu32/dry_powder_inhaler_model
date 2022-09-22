/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_ELECTROSTATIC_SCREEN,electrostatic/screened,14)
#else

#ifndef COHESION_ELECTROSTATIC_SCREENED_H_
#define COHESION_ELECTROSTATIC_SCREENED_H_

#include "contact_models.h"
#include "math.h"
#include "math_extra_liggghts.h"
#include "global_properties.h"
#include "fix_charge_gran.h"
#include "fix_property_particle.h"

namespace LIGGGHTS {

namespace ContactModels {
   

  template<>
  class CohesionModel<COHESION_ELECTROSTATIC_SCREEN> : protected Pointers 
  {
  
  protected :
  
  double permittivity;
  double screeningDistance;
  int history_offset;
  FixPropertyParticle * fix_charge;
  FixPropertyParticle * fix_polarization;
  bool longRangeFlag;
  
  #include "tensorOperations.h"
  
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup) :
      Pointers(lmp), permittivity(0.0), screeningDistance( 0.0 )
    {
        fprintf( screen, "Screened electrostatic force model active.\n" );  
	fix_charge = NULL;
        fix_polarization = NULL;
        longRangeFlag = true;
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) 
    {
    
       // -- vacuum permittivity --
       registry.registerProperty("permittivity", &MODEL_PARAMS::createPermittivity);      
       registry.connect("permittivity", permittivity,"cohesion_model_electrostatic_screened");
       
       // -- screening distance for P3M --
       registry.registerProperty("screeningDistance", &MODEL_PARAMS::createScreeningDistance);      
       registry.connect("screeningDistance", screeningDistance,"cohesion_model_electrostatic_screened");
	
       //int i = modify->find_fix("efieldtransfer");
       //if (i < 0) error->all(FLERR,"Illegal efieldtransfer command, need a fix called 'efieldtransfer'");            
       //fix_FixEfieldGran_ = static_cast<FixEfieldGran*>(modify->fix[i]);
	
       void* i = modify->find_fix_property("charge","property/particle","scalar",0,0,"efield/gran");
      
       if( !i ) error->all(FLERR,"Atoms needs to have property 'charge'!");
       fix_charge = static_cast<FixPropertyParticle*>(modify->find_fix_property("charge","property/particle","scalar",0,0,"efield/gran"));

       fix_polarization = static_cast<FixPropertyParticle*>(modify->find_fix_property("polarization","property/particle","vector",0,0,"efield/polarization",false));
       
       if( fix_polarization ) fprintf( screen, "Dielectrophoresis active...\n" );
         
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

    // compute the electric field due to charge on particle and gradient of electric field due to 
    void computeEfield( const double r, int i, const Tensor& n, Tensor& efield, Tensor& gradEfield ) const
    {
	
	// Charge on the particle
	const double charge_i = (fix_charge->vector_atom)[i];
        
	const double pi = M_PI;
		
	
	const double r2 = r*r;
	const double sqrt2sigma = sqrt( 2.0 ) * screeningDistance;
	const double sigma2 = 2.0 * screeningDistance * screeningDistance;
	
	const double erfcExp = erfc( r/sqrt2sigma );
	const double gaussExp = exp( -r2/sigma2 );
	
		
        const double Fn_coh = 1./(4. * pi * permittivity * r) *
		 		( 
				    erfcExp / r 
				  + sqrt( 2.0/pi ) * gaussExp / screeningDistance 
				 ); 
	
	const double gradFn_coh = 1./( 4.*pi*permittivity ) * 
				 ( 
				     sqrt( 2.0/pi ) * gaussExp / ( screeningDistance * r * r ) 
				   + sqrt( 2.0/pi ) * gaussExp / ( screeningDistance * screeningDistance * screeningDistance )
				   + erfcExp / ( r2 * r )
				 ); 
	
	// -- add Coulombic field contribution to electric field --
	efield += n * ( Fn_coh * charge_i );
	
	const Tensor n2 = n * ( n.transpose() );
	
	// -- add Coulombic field contribution to electric field gradient -- 
	gradEfield -= ( ( Tensor::eye() - n2 ) * ( Fn_coh * charge_i ) ); 
	gradEfield -= ( n2 * ( gradFn_coh * charge_i ) );        
	
	// -- Polarization Contribution --
	if( !fix_polarization ) return; // if no electric dipoles are defined, return
	
	
	// -- electric dipole of particle i --
	const Tensor pii( fix_polarization->array_atom[i] );
	
	// -- update dipole contribution to the electric field --
	efield += ( ( Tensor::eye() - n2 ) * pii * Fn_coh );
	efield -= ( n2 * pii * gradFn_coh ); 					   

	const double gradFp_coh1 = 1.0/(4.*pi*permittivity) * 
				   (
				       3 * erfcExp/ (r2*r2)
				     + sqrt( 2.0/pi ) * gaussExp / ( screeningDistance * screeningDistance * screeningDistance * r )
				     + 3.0 * sqrt( 2.0/pi ) * gaussExp / ( screeningDistance * r2 * r )
				   );
	
	
	const double gradFp_coh2 = 1.0/(4.*pi*permittivity) * 
		 		   (
				       9.0 * erfcExp/ (r2*r2)
				     + sqrt( 2.0/pi ) * gaussExp / ( screeningDistance * screeningDistance * screeningDistance * r )
				     + 9.0 * sqrt( 2.0/pi ) * gaussExp / ( screeningDistance * r2 * r )
				   );
	
	const double gradFp_coh3 = 1.0/(2.*pi*permittivity) * 
				   (
				       sqrt( 2.0/pi ) * gaussExp * r / ( pow(screeningDistance,5) )
				   );
	
	
	const double pin = ( pii.transpose() * n ).value();
	
	const Tensor Psd1 = n2 * pin;

	const Tensor Psd2 =   Tensor::eye() * pin 
		      	    + pii * ( n.transpose() ) 
		            - Psd1 * 2.0;  
	
	
	gradEfield -= (  Psd2 * gradFp_coh1
	 	      +	( pii * ( n.transpose() ) ) * gradFp_coh1
		      - Psd1 * ( gradFp_coh2 + gradFp_coh3 ) ); 
	 
    }
    

    void computeEfieldCollision( CollisionData & cdata, Tensor& efield, Tensor& gradEfield ) const
    {
	
	// -- normal vector --
	const Tensor n( cdata.en );
	const double r = sqrt( cdata.rsq );
	const int i = cdata.i;
        
	computeEfield( r, i, n, efield, gradEfield );
	
    }
    
    void computeEfieldnoCollision( ContactData & cdata, Tensor& efield, Tensor& gradEfield ) const
    {
        
	const double r = sqrt( cdata.rsq );
	const int i = cdata.i;
	
	Tensor n( cdata.delta, false );
	n *= (1.0/r);
	
	computeEfield( r, i, n, efield, gradEfield );
	
    }
    
    
    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
              
	      
       
       // -- electrif field at the particle j position due to particle i --
       Tensor efield( 3, 1 );
       
       // -- electrif field gradient at the particle j position due to particle i --
       Tensor gradEfield( 3, 3 );
       
       // -- compute electric field at particle j position --
       computeEfieldCollision( cdata, efield, gradEfield ); 
       
       const int i = cdata.i;
       const double charge_i = (fix_charge->vector_atom)[i];
       int j = cdata.is_wall ? cdata.i : cdata.j;
       	
       double charge_j = cdata.is_wall ? -charge_i : (fix_charge->vector_atom)[j];	
	
	// -- force acting on particle j due to particle i --
       Tensor fn = efield * charge_j;
       
       if( fix_polarization )
       { 
	   Tensor pj( fix_polarization->array_atom[j] ); // if wall, j = i -> pj = pi (mirroring)
	   fn -= ( gradEfield * pj );   
	   	   
	   //fprintf( screen, "Collision: fn = (%e %e %e) \n ", fn[0], fn[1], fn[2] ); 
       }	
	
       if( !cdata.is_wall )
       {
           for( int k = 0; k < 3; ++k )
	   {
              i_forces.delta_F[k] += fn[k];
              j_forces.delta_F[k] -= fn[k];  
	   }
       }else{
           for( int k = 0; k < 3; ++k )
	   {
              i_forces.delta_F[k] += fn[k];
	   }   
       }	
	 
       cdata.has_force_update = true || cdata.has_force_update;
	      
    }

    void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces) 
    {
       
       // -- electrif field at the particle j position due to particle i --
       Tensor efield( 3, 1 );
       
       // -- electrif field gradient at the particle j position due to particle i --
       Tensor gradEfield( 3, 3 );
       
       // -- compute electric field at particle j position --
       computeEfieldnoCollision( cdata, efield, gradEfield ); 
       
       const int i = cdata.i;
       const double charge_i = (fix_charge->vector_atom)[i];
       int j = cdata.is_wall ? cdata.i : cdata.j;
       	
       double charge_j = cdata.is_wall ? -charge_i : (fix_charge->vector_atom)[j];	
	
	// -- force acting on particle j due to particle i --
       Tensor fn = efield * charge_j;
       
       if( fix_polarization )
       { 
	   Tensor pj( fix_polarization->array_atom[j] ); // if wall, j = i -> pj = pi (mirroring)
	   fn -= gradEfield * pj;  
	   
	   //fprintf( screen, "NoCollision: fn = (%e %e %e) \n ", fn[0], fn[1], fn[2] );  
       }	
	
       if( !cdata.is_wall )
       {
           for( int k = 0; k < 3; ++k )
	   {
              i_forces.delta_F[k] += fn[k];
              j_forces.delta_F[k] -= fn[k];  
	   }
       }else{
           for( int k = 0; k < 3; ++k )
	   {
              i_forces.delta_F[k] += fn[k];
	   }   
       }	
	 
       cdata.has_force_update = true || cdata.has_force_update;
       
    }
      
  }; // -- end of cohesion class --
    
};

}

#endif // COHESION_ELECTROSTATIC_SCREENED_H_
#endif


















