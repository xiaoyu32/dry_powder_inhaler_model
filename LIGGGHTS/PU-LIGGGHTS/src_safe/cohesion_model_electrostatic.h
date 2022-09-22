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
COHESION_MODEL(COHESION_ELECTROSTATIC,electrostatic,3)
#else

#ifndef COHESION_MODEL_ELECTROSTATIC_H_
#define COHESION_MODEL_ELECTROSTATIC_H_

#include "contact_models.h"
#include "math.h"
#include "math_extra_liggghts.h"
#include "global_properties.h"
#include "fix_charge_gran.h"
#include "fix_property_particle.h"

namespace LIGGGHTS {

namespace ContactModels {

  template<>
  class CohesionModel<COHESION_ELECTROSTATIC> : protected Pointers {
  
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup) :
      Pointers(lmp), permittivity(0.0)
    {
    
        fix_polarization = NULL;
        longRangeFlag = true;
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) 
    {
       registry.registerProperty("permittivity", &MODEL_PARAMS::createPermittivity);      
       registry.connect("permittivity", permittivity,"cohesion_model_electrostatic");
	
       //int i = modify->find_fix("efieldtransfer");
       //if (i < 0) error->all(FLERR,"Illegal efieldtransfer command, need a fix called 'efieldtransfer'");            
       //fix_FixEfieldGran_ = static_cast<FixEfieldGran*>(modify->fix[i]);
	
       void* i = modify->find_fix_property("charge","property/particle","scalar",0,0,"efield/gran");
      
       if( !i ) error->all(FLERR,"Atoms needs to have property 'charge'!");
       fix_charge = static_cast<FixPropertyParticle*>(modify->find_fix_property("charge","property/particle","scalar",0,0,"efield/gran"));

       fix_polarization = static_cast<FixPropertyParticle*>(modify->find_fix_property("polarization","property/particle","vector",0,0,"efield/polarization",false));
       fix_cell_center = static_cast<FixPropertyParticle*>(modify->find_fix_property("cell_center","property/particle","vector",3,0,"couple/cfd/electric",false));
       
       
       void* ptr = modify->find_fix_property("electricfield_gradient","property/particle","vector",9,0,"couple/cfd/electric", false); 
       
       if( ptr ){
           longRangeFlag = true;
       }else{
           longRangeFlag = false;
       }
	
       
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    
    void cross_product( double* a, double* b, double* c )
    {
        
	c[0] = a[1]*b[2] - b[1]*a[2];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
	
    }
    
    void dielectrophoresis_collision( CollisionData & cdata, ForceData & i_forces, ForceData & j_forces, int i, int j, int sign )
    {
        
    
        double pi[3];
	double pj[3];
	double dd[3];
        
	double r = sqrt( cdata.rsq );
	
	const double charge_i = (fix_charge->vector_atom)[i];
        
	double charge_j;
	
	if( !cdata.is_wall ){
	    
	    charge_j = (fix_charge->vector_atom)[j];
	    
	    for( int k = 0; k < 3; ++k )
	    { 
	        pj[k] = fix_polarization->array_atom[j][k];
		dd[k] = sign*cdata.en[k]*r;
	    }
	    
	}else{
	
	    charge_j = -charge_i;
	    
	    for( int k = 0; k < 3; ++k )
	    { 
	        pj[k] = fix_polarization->array_atom[i][k];
		dd[k] = sign*cdata.en[k]*r;
	    }
	    
	}
	
	
	for( int k = 0; k < 3; ++k )
	{ 
	    pi[k] = fix_polarization->array_atom[i][k];
	}

	double dE[3][3];
        double E[3];
	double qE[3];
	
	compute_grade( E, qE, dE, pj, charge_j, dd );
        
	double ff[3];

	for( int k = 0; k < 3; ++k )
	   ff[k] = 0.0;

	for( int k1 = 0; k1 <3; ++k1 )
	{
	    for( int k2 = 0; k2 < 3; ++k2 )
	    {
		//gradient of electric field
		ff[k1] += dE[k1][k2]*pi[k1];
		  	
	    }
	    
	    
	    //electric field contribution
	    ff[k1] +=  charge_i * E[k1];
	    
	    E[k1] += qE[k1];
	    
	}
	
	double tor[3];
	
	cross_product( pi, E, tor );
	
	/*
	fprintf( screen, "\n" );
	fprintf( screen, "collided! \n" );
	fprintf( screen, "dd = %e %e %e \n", dd[0], dd[1], dd[2] );
	fprintf( screen, "pi = %e %e %e \n", pi[0], pi[1], pi[2] );
	fprintf( screen, "pj = %e %e %e \n", pj[0], pj[1], pj[2] );
	fprintf( screen, "E = %e %e %e \n", E[0], E[1], E[2] );
	fprintf( screen, "no T = %e %e %e \n", tor[0], tor[1], tor[2] );
        */
	/*fprintf( screen, "\n" );
	fprintf( screen, "id = %d \n", i );
        fprintf( screen, "x = %e %e %e \n", atom->x[i][0], atom->x[i][1], atom->x[i][2] );
	fprintf( screen, "F = %e %e %e \n", ff[0], ff[1], ff[2] );
	fprintf( screen, "\n" );*/
	
	
        
	if( !cdata.is_wall ){
	
	     i_forces.delta_F[0] += ff[0];
      	     i_forces.delta_F[1] += ff[1];
      	     i_forces.delta_F[2] += ff[2];
	     
	     j_forces.delta_F[0] -= ff[0];
      	     j_forces.delta_F[1] -= ff[1];
      	     j_forces.delta_F[2] -= ff[2];   
	     
	     i_forces.delta_torque[0] += tor[0];
             i_forces.delta_torque[1] += tor[1];
             i_forces.delta_torque[2] += tor[2];
	     
	     j_forces.delta_torque[0] -= tor[0];
             j_forces.delta_torque[1] -= tor[1];
             j_forces.delta_torque[2] -= tor[2];
	     
	}else{
	
	     i_forces.delta_F[0] += ff[0];
      	     i_forces.delta_F[1] += ff[1];
      	     i_forces.delta_F[2] += ff[2];
	     
	     
	     i_forces.delta_torque[0] += tor[0];
             i_forces.delta_torque[1] += tor[1];
             i_forces.delta_torque[2] += tor[2];
	     	
	}
	
    }
    
    void dielectrophoresis( ContactData & cdata, ForceData & i_forces, ForceData & j_forces, int i, int j, int sign )
    {
        
    	
        double pi[3];
	double pj[3];
	double dd[3];
        
	double r = sqrt( cdata.rsq );
	
	const double charge_i = (fix_charge->vector_atom)[i];
        
	double charge_j;
	
	if( !cdata.is_wall ){
	    
	    	    
	    charge_j = (fix_charge->vector_atom)[j];
	    
	    for( int k = 0; k < 3; ++k )
	    { 
	        pj[k] = fix_polarization->array_atom[j][k];
		dd[k] = atom->x[i][k] - atom->x[j][k];
	    }
	    
	}else{
	
	    charge_j = -charge_i;
	    
	    for( int k = 0; k < 3; ++k )
	    { 
	        pj[k] = fix_polarization->array_atom[i][k];
		dd[k] = sign*cdata.delta[k];
	    }
	    
	}
	
	
	for( int k = 0; k < 3; ++k )
	{ 
	    pi[k] = fix_polarization->array_atom[i][k];
	}

	double dE[3][3];
        double E[3];
	double qE[3];
	
	compute_grade( E, qE, dE, pj, charge_j, dd );
        
	double ff[3];

	for( int k = 0; k < 3; ++k )
	   ff[k] = 0.0;

	for( int k1 = 0; k1 <3; ++k1 )
	{
	    for( int k2 = 0; k2 < 3; ++k2 )
	    {
		//gradient of electric field
		ff[k1] += dE[k1][k2]*pi[k1];
		  	
	    }
	    
	    
	    //electric field contribution
	    ff[k1] +=  charge_i * E[k1];
	    
	    E[k1] += qE[k1];
	    
	}
	
	double tor[3];
	
	cross_product( pi, E, tor );
	
	/*
	fprintf( screen, "\n" );
	fprintf( screen, "dd = %e %e %e \n", dd[0], dd[1], dd[2] );
	fprintf( screen, "pi = %e %e %e \n", pi[0], pi[1], pi[2] );
	fprintf( screen, "pj = %e %e %e \n", pj[0], pj[1], pj[2] );
	fprintf( screen, "E = %e %e %e \n", E[0], E[1], E[2] );
	fprintf( screen, "no T = %e %e %e \n", tor[0], tor[1], tor[2] );
        fprintf( screen, "F = %e %e %e \n", ff[0], ff[1], ff[2] );
	fprintf( screen, "\n" );
	*/
	
	/*fprintf( screen, "\n" );
	fprintf( screen, "id = %d \n", i );
	fprintf( screen, "x = %e %e %e \n", atom->x[i][0], atom->x[i][1], atom->x[i][2] );
	fprintf( screen, "F = %e %e %e \n", ff[0], ff[1], ff[2] );
	fprintf( screen, "\n" );*/
	
	
	if( !cdata.is_wall ){
	
	     i_forces.delta_F[0] += ff[0];
      	     i_forces.delta_F[1] += ff[1];
      	     i_forces.delta_F[2] += ff[2];
	     
	     j_forces.delta_F[0] -= ff[0];
      	     j_forces.delta_F[1] -= ff[1];
      	     j_forces.delta_F[2] -= ff[2];   
	     
	     
	     i_forces.delta_torque[0] += tor[0];
             i_forces.delta_torque[1] += tor[1];
             i_forces.delta_torque[2] += tor[2];
	     
	     j_forces.delta_torque[0] -= tor[0];
             j_forces.delta_torque[1] -= tor[1];
             j_forces.delta_torque[2] -= tor[2];
	     
	     
	}else{
	
	     i_forces.delta_F[0] += ff[0];
      	     i_forces.delta_F[1] += ff[1];
      	     i_forces.delta_F[2] += ff[2];
	     
	     
	     i_forces.delta_torque[0] += tor[0];
             i_forces.delta_torque[1] += tor[1];
             i_forces.delta_torque[2] += tor[2];
	     
	     
	}
	
    }
    
    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
              
       //dielectrophoresis
       if( fix_polarization )
       {
           int i = cdata.i;
	   int j;
	   
	   if( !cdata.is_wall )
	   {
	      j = cdata.j;
	   }else{
	      j = i;
	   }
       
	   dielectrophoresis_collision( cdata, i_forces, j_forces, i, j, 1 );
	   //if( !cdata.is_wall ) dielectrophoresis_collision( cdata, j_forces, i_forces, j, i, -1 );
	   
	   
       }
        	
       if( !cdata.is_wall ){
         // The distance between the sphere's centeres
         const double r = sqrt( cdata.rsq );
         const int i = cdata.i;
         const int j = cdata.j;      
         const double pi = M_PI;
      
         //if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
         //double * const shear = &cdata.contact_history[history_offset];
   
         const double charge_i = (fix_charge->vector_atom)[i];
         const double charge_j = (fix_charge->vector_atom)[j];
	              
         double Fn_coh;
         Fn_coh = 1./(4.*pi*permittivity)*charge_i*charge_j/(r*r); 	
         //cdata.Fn += Fn_coh;
	 
	 //fprintf( logfile, "Fn_coh(2) = %e \n", Fn_coh );
	 
         // apply normal force
         const double fx = Fn_coh * cdata.en[0];
         const double fy = Fn_coh * cdata.en[1];
         const double fz = Fn_coh * cdata.en[2];

         i_forces.delta_F[0] += fx;
         i_forces.delta_F[1] += fy;
         i_forces.delta_F[2] += fz;

         j_forces.delta_F[0] -= fx;
         j_forces.delta_F[1] -= fy;
         j_forces.delta_F[2] -= fz;

         // store for noCollision
         //shear[0] = 1.0;
       }else{
       
         //conductive walls
	  
	 double r = sqrt( cdata.rsq );
	 const int i = cdata.i;
         const double pi = M_PI;
      
         //if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
         //double * const shear = &cdata.contact_history[history_offset];
   
         const double charge_i = (fix_charge->vector_atom)[i];
         double Fn_coh;
         Fn_coh = -1./(4.*pi*permittivity)*charge_i*charge_i/(4*r*r);
	 
	 //fprintf( screen, "Charge1: %e \n", Fn_coh );
	  	
         //cdata.Fn += Fn_coh;

         // apply normal force
         const double fx = Fn_coh * cdata.en[0];
         const double fy = Fn_coh * cdata.en[1];
         const double fz = Fn_coh * cdata.en[2];

         i_forces.delta_F[0] += fx;
         i_forces.delta_F[1] += fy;
         i_forces.delta_F[2] += fz;
	 
	 //correction to avoid double counting the wall	 
	 if( fix_cell_center )
	 {
	     
	    /*double* cell =  fix_cell_center->array_atom[i];
            
	    double wall_delta = 0;
	    double wallpos[3];
	    double rdelta = 0;
	    
	    for( int k = 0; k < 3; ++k )
	    {
	    	rdelta += cdata.delta[k] * cdata.delta[k];
	    }
	    
	    rdelta = sqrt( rdelta );
	    
	    for( int k = 0; k < 3; ++k )
	    {
	        wallpos[k] = cdata.delta[k]*r/rdelta + atom->x[i][k] -cell[k];
	    }
	    
	    for( int k = 0; k < 3; ++k )
	    {
	    	wall_delta += cdata.delta[k]/rdelta * wallpos[k];
	    }
	    	    
	    double efield = charge_i/(wall_delta*wall_delta)/(4.*pi*permittivity);
	     
	    const double Fn_cor = charge_i*efield;
	    
	    // apply normal force
            const double fxc = Fn_cor * cdata.en[0];
            const double fyc = Fn_cor * cdata.en[1];
            const double fzc = Fn_cor * cdata.en[2];
	 
	    i_forces.delta_F[0] += fxc;
            i_forces.delta_F[1] += fyc;
            i_forces.delta_F[2] += fzc;*/
	     
	 }
	 
	 
       }
       
       cdata.has_force_update = true||cdata.has_force_update;
       
    }

    void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces) 
    {
        
	//fprintf( screen, "No Collision \n" );
	      
       //dielectrophoresis
       if( fix_polarization  )
       {
	   
	   int i = cdata.i;
	   int j;
	   
	   if( !cdata.is_wall )
	   {
	      j = cdata.j;
	   }else{
	      j = i;
	   }
	   
	   dielectrophoresis( cdata, i_forces, j_forces, i, j, 1 );
	   	   
	   if( !cdata.is_wall ) 
	   {
	   //    dielectrophoresis( cdata, j_forces, i_forces, j, i, 1 );
	   }  
       }
       
       if( !cdata.is_wall ){ 
         
	 
	 double dd[3];   
	 const double r = sqrt(cdata.rsq);
         const double rinv = 1/r;
         double rr = 0;
	 const int i = cdata.i;
         const int j = cdata.j;	  
         const double pi = M_PI;
      
         //if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
         //double * const shear = &cdata.contact_history[history_offset];

         const double charge_i = (fix_charge->vector_atom)[i];
         const double charge_j = (fix_charge->vector_atom)[j];
	 
	 for( int k = 0; k < 3; ++k )
	 {
	     dd[k] = atom->x[j][k] - atom->x[i][k];
	     rr += dd[k]*dd[k];
	 } 
	 rr = 1.0/sqrt(rr); 
	              
         double Fn_coh;
         Fn_coh = -1./(4.*pi*permittivity) * charge_i* charge_j * rr * rr;
	 //printf( "Fn_coh(1) = %e \n", Fn_coh );
	 
         const double fx = Fn_coh * dd[0] * rr;
         const double fy = Fn_coh * dd[1] * rr;
         const double fz = Fn_coh * dd[2] * rr;
	 
   	 i_forces.delta_F[0] += fx;
      	 i_forces.delta_F[1] += fy;
      	 i_forces.delta_F[2] += fz;

      	 j_forces.delta_F[0] -= fx;
      	 j_forces.delta_F[1] -= fy;
      	 j_forces.delta_F[2] -= fz;      
        
	
       }else{
         
	 //FIXME: never active with granular walls
	 
	 //conductive walls
	 double r = sqrt( cdata.rsq );
	 const double rinv = 1.0/r;
	 const int i = cdata.i;
         const double pi = M_PI;
      
         const double charge_i = (fix_charge->vector_atom)[i];
         double Fn_coh;
         Fn_coh = -1./(4.*pi*permittivity)*charge_i*charge_i/(4*r*r);
	 
	 //fprintf( screen, "hoppsan: %e \n", Fn_coh );
	  	
         //cdata.Fn += Fn_coh;

         // apply normal force
         const double fx = Fn_coh * rinv * cdata.delta[0];// * cdata.en[0];
         const double fy = Fn_coh * rinv * cdata.delta[1];// * cdata.en[1];
         const double fz = Fn_coh * rinv * cdata.delta[2];// * cdata.en[2];

         i_forces.delta_F[0] += fx;
         i_forces.delta_F[1] += fy;
         i_forces.delta_F[2] += fz;
	 
       }
       
       cdata.has_force_update = true||cdata.has_force_update;
    }

  private:
    
    
    //compute gradient of electric field for dielectrophoresis
    void compute_grade( double* E, double* qE, double dE[][3], double* p, double q, double* x )
    {
       
	double r = 0;
	double px = 0;
	
	double cof = 1.0/(4.0*M_PI*permittivity);
	
	for( int i = 0; i < 3; ++i )
	{
	    r += x[i]*x[i];
	    px += x[i]*p[i];
	}
	
	r = sqrt(r);
	
	if( r <= 0 ) return;
		
	for( int i = 0; i < 3; ++i )
	{
	
	   //fprintf( screen, "px=%e \n", px );
	   E[i] = cof*( 3*px*x[i]/(r*r*r*r*r) - p[i]/(r*r*r));
	   
	   qE[i] = -cof*q*x[i]/(r*r*r);
	   
	   
	   for( int j = 0; j < 3; ++j )
	   {
	       
	       dE[i][j] =  3*x[i]*p[j]/(r*r*r*r*r)
			  +3*p[i]*x[j]/(r*r*r*r*r)
			  -15*px*x[i]*x[j]/(r*r*r*r*r*r*r);

	       if( i == j )
	       {
		   dE[i][j] += 3*px/(r*r*r*r*r);
	       }

		
	       //charge contribution to the electric gradient (only take into account if no long range information is available	
	       if( !longRangeFlag || true )
	       {
		

	       	    dE[i][j] += -3*q/(r*r*r*r*r)*x[i]*x[j];

		    if( i == j )
		    {
			dE[i][j] += q/(r*r*r);
		    }
	       	  
	       }
	       		  
	       dE[i][j] *= cof;
	       
	   }
	}
	
    }
        
    double permittivity;
    int history_offset;
    FixPropertyParticle* fix_charge;
    FixPropertyParticle* fix_polarization;
    FixPropertyParticle* fix_cell_center;
    bool longRangeFlag;
    
    
  };
}

}

#endif // COHESION_MODEL_ELECTROSTATIC_H_
#endif
