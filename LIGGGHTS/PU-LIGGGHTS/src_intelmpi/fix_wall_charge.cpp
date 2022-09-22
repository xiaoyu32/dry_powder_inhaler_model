#include "math.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "mech_param_gran.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair_gran.h"

#include "fix_wall_charge.h"
#include "fix_charge_gran.h"

#include "efield_model.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixWallCharge::FixWallCharge(class LAMMPS *lmp, int narg, char **arg) : FixWallGran(lmp, narg, arg)
{
  

  int iarg = 5;

  area_correction_flag = 0;
  
  transfer_acceleration = 1.0;
  
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg],"area_correction") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'area_correction'");
      if(strcmp(arg[iarg+1],"yes") == 0)
        area_correction_flag = 1;
      else if(strcmp(arg[iarg+1],"no") == 0)
        area_correction_flag = 0;
      else error->fix_error(FLERR,this,"");
      iarg += 1;
      hasargs = true;
    } 

    ++iarg;
    
  }
  
  fix_efieldgran = NULL;
  fix_work_wall = NULL;
  work_wall = NULL;
  fix_p = NULL;
  fix_resistivity_wall = NULL;
  resistivity_wall = NULL;
}


void FixWallCharge::post_create(){
    
    char * fixarg[3];
    
    fix_efieldgran = static_cast<FixEfieldGran*> ( modify->find_fix_id("wall_gran_efield") );
    
    if( !fix_efieldgran ){
       fixarg[0] = (char*) "wall_gran_efield";
       fixarg[1] = (char*) "all";
       fixarg[2] = (char*) "efield/gran";
    
       modify->add_fix(3,fixarg);
    }
  
    FixWallGran::post_create();
    
}

FixWallCharge::~FixWallCharge()
{
  if( work_wall ) delete [] work_wall;
  if( resistivity_wall ) delete [] resistivity_wall;
}


void FixWallCharge::init()
{

   fix_efieldgran = static_cast<FixEfieldGran*> ( modify->find_fix_id("wall_gran_efield") );
   fix_work_wall = static_cast<FixPropertyGlobal*>(modify->find_fix_property("work_wall","property/global","scalar",0,0,style));
   
   if( !fix_work_wall )
      error->fix_error(FLERR,this,"Fix wall/gran/efield requires work function value for walls!");

   if(work_wall) delete [] work_wall;
   //TODO: implement more wall types
   work_wall = new double[1];
   work_wall[0] = fix_work_wall->compute_scalar();
   
   if (fix_efieldgran->get_relaxation_model_flag())
   {   
   	fix_resistivity_wall = static_cast<FixPropertyGlobal*>(modify->find_fix_property("resistivity_wall","property/global","scalar",0,0,style));

   	if( !fix_resistivity_wall )
      		error->fix_error(FLERR,this,"Fix wall/gran/efield requires resistivity  value for walls!");

   	if(resistivity_wall) delete [] resistivity_wall;
        resistivity_wall = new double[1];
        resistivity_wall[0] = fix_resistivity_wall->compute_scalar();
   }

       
   FixWallGran::init();
   
   fix_p = static_cast<FixPropertyParticle*>(modify->find_fix_property("permanent_polarization","property/particle","vector",3,0,this->style,false));
   
}

void FixWallCharge::post_force(int vflag)
{
   FixWallGran::post_force(vflag);
}


int FixWallCharge::setmask()
{
  int mask = FixWallGran::setmask();
  mask |= POST_FORCE;
  return mask;
}

double FixWallCharge::abs( double a )
{
   if( a > 0 )  return a;
   else		return -a;
}

//Overwrite collision model
void FixWallCharge::compute_force(LCM::CollisionData & cdata, double *vwall)
{
   
   
   double cA, cA_old;
   double flux, dirFlux[3], flux_old;
   double reff, electricf, electricf_ext;
   double q, dA;
   double wa, wb, dwork;
   double w_cof_a, w_cof_b;
   double delta_old[3];
   double delta[3];
   double deltan_old;
   double foo,bar;
   double r;
   int iPart = cdata.i;
   
   int iType = atom->type[iPart];
   int *mask = atom->mask;
   double idt = 1.0/update->dt;   
   
   FixPropertyGlobal* fix_acc = static_cast<FixPropertyGlobal*>(modify->find_fix_property("acceleration","property/global","scalar",0,0,style,false));
   
   if( fix_acc ) transfer_acceleration = fix_acc->compute_scalar();
   
   const double chiStiff = 1.0/lmp->force->chiStiffnessScaling(); 

   if( !( mask[iPart] & groupbit )  ) return;

   //double ri = cdata.radi;
   //double rj = cdata.radj;
   double delx, dely, delz;
   int charging_model_type;
   
   fix_efieldgran->updatePtrs();

   double perm = fix_efieldgran->permittivity;
   double delta_charge = fix_efieldgran->delta_charge;
   double bd_field = fix_efieldgran->bd_field;
	   
   // normal vector
   delx = cdata.delta[0];
   dely = cdata.delta[1];
   delz = cdata.delta[2];
   
   if( !cdata.is_wall ){
      return;
   }  

   // no collision TODO: implement cut-off distance
   if( cdata.deltan <= 0 ) return;
   
   //reference radius
   //reff = cdata.is_wall ? ri : (ri*rj/(ri+rj));
   reff = cdata.r;  
   r = atom->radius[iPart];
   //particles current charge
   q = fix_efieldgran->charge[iPart];
   
   foo = 0.0;
   //normal vector pointing from the contact point to particle center
   
   for( int k = 0; k < 3; ++k )
   {
       foo += cdata.delta[k]*cdata.delta[k];
   }
   
   foo = sqrt(foo);
   
   double normalVector[3];
   
   for( int k = 0; k < 3; ++k )
       normalVector[k] = cdata.delta[k]/foo;
   
   electricf_ext = 0.0;
   
   //external and ambient electricfield contribution
   if( fix_efieldgran->ef_coupling_flag != 0 )
   {
	    
       //linear interpolation for intersection point electricfield
       for( int k = 0; k <3; ++k )
       {
          electricf_ext += cdata.delta[k]/foo * fix_efieldgran->fix_ef_coupling->array_atom[iPart][k];
       }
	
   }
   
   
   
   //old contact area 
   if( fix_efieldgran->prev_loc_flag ){

      for( int k = 0; k < 3; ++k ){
	  delta[k] = cdata.delta[k]*(r-cdata.deltan)/foo;
      }

      deltan_old = 0.0;  
      for( int k = 0; k < 3; ++k )
      {
          delta_old[k] = delta[k] + (fix_efieldgran->prev_loc[iPart][k] - atom->x[iPart][k]);
      }

      //project to wall normal to account change of contact point
      deltan_old = 0.0;
      
      for( int k = 0; k < 3; ++k )
      {
           deltan_old += delta_old[k]*delta[k];
      }
      
      deltan_old /= (r-cdata.deltan);

      deltan_old = r-deltan_old;
      cA_old = M_PI * r * deltan_old;

      //no collision on the previous time step 	
      if( cA_old < 0.0 )
      { 
         cA_old = 0.0;
      }

   }else{
      cA_old = 0.0;
   }   

   
   //current contact area 
   cA = M_PI * r * cdata.deltan;

   dA = cA-cA_old;
   w_cof_a = fix_efieldgran->work_a[iType-1];
   w_cof_b = fix_efieldgran->work_b[iType-1];
   charging_model_type = fix_efieldgran->work_model[iType-1];   

   wa = FixEfieldGran::work( reff, w_cof_a, w_cof_b, charging_model_type );
   wb = FixEfieldGran::work( 1.0 , work_wall[0], 0, WALL );

   dwork = perm / (delta_charge*ELECTRON_CHARGE ) * (wa - wb);
  
   /* 
   if( dA < 0 )
   {
      dA = 0.0;
      return;   
   }
   */
   
   
   flux = 0.0;
   flux_old = 1.0;
   
   int it = 0;
   
   //implicit Euler for flux
   //electricf = -2*q/(r*r);
   //electricf /= (4*M_PI*perm);

   double polField = 0;

   if( fix_p )
   {
      double polarizationEfield[3] = {0,0,0};
      fix_efieldgran->efieldModel()->computeEfieldPolarization( normalVector, r, fix_p->array_atom[iPart], polarizationEfield );

      for( int k = 0; k < 3; ++k ) // -- mirroring --
          polField += 2.0 * polarizationEfield[k] * normalVector[k];
   }
   
   electricf = fix_efieldgran->efieldModel()->implElectricField( q, -q, r, r );
   
   //double impl_coff = 1.0 + 2 * chiStiff * transfer_acceleration * dA/(4*M_PI*r*r);
   double impl_coff = 1.0 + chiStiff * transfer_acceleration * dA * fix_efieldgran->efieldModel()->implCoeff( r, r ); 
   
   if (dA > 0)
   {
	flux = chiStiff * transfer_acceleration * idt * dA *( dwork - perm * ( electricf + electricf_ext + polField ) );
	if (fix_efieldgran->get_relaxation_model_flag())
	{
		double part_resistivity = fix_efieldgran->resistivity[iType-1];
		double wall_resistivity = resistivity_wall[0];
		double effective_resistivity = part_resistivity > wall_resistivity? part_resistivity: wall_resistivity;
		flux -= chiStiff*transfer_acceleration*idt*(electricf + electricf_ext)/effective_resistivity*cA;		
	}
   }
   else if (dA < 0)
   {
	flux = 0.0;
	if (fix_efieldgran->get_relaxation_model_flag())
        {
                double part_resistivity = fix_efieldgran->resistivity[iType-1];
                double wall_resistivity = resistivity_wall[0];
                double effective_resistivity = part_resistivity > wall_resistivity? part_resistivity: wall_resistivity;
                flux -= chiStiff*transfer_acceleration*idt*(electricf + electricf_ext)/effective_resistivity*cA;
        }
   }
   else
   {
	flux = 0.0;
   }
	
   flux /= impl_coff; 
   
   
   //conductive walls by mirroring, assume deltan << reff
   //electricf = -2*(q-flux/idt)/(r*r);
   //electricf /= (4*M_PI*perm); 
   
   electricf = fix_efieldgran->efieldModel()->implElectricField( q-flux/idt, -q+flux/idt, r, r );
   
   if (  fix_efieldgran->bd_model_flag()  )
   {
   	//fprintf( screen, "Wall is using new breakdown model...");
	
	if( abs( electricf + electricf_ext + polField ) > bd_field )
	{
		double surface_area = 4.0*M_PI*r*r;
       
           	if ( (electricf + electricf_ext) >= 0 )
           	{
       	        	double q_12 = (bd_field - (electricf + electricf_ext + polField))/(2.0/surface_area)*perm;
       	        	//flux = 1.1*q_12*idt;
			flux = q_12*idt;
           	}
       
           	else  
           	{
                	double q_12 = -(bd_field + (electricf + electricf_ext + polField))/(2.0/surface_area)*perm;
       	        	//flux = 1.1*q_12*idt; 
			flux = q_12*idt;     
           	}
       
       
       		fprintf( screen, "New break down at the wall! (electricf= %e) (charge = %e)\n", electricf+electricf_ext, q  );
       		fprintf( screen, "Work = %e electric_local = %e electric = %e \n", dwork, perm*electricf, perm*electricf_ext );
	}
	  
   }
   else
   {
   	//fprintf( screen, "Wall is using old breakdown model...");
	
	if( abs( electricf + electricf_ext ) > bd_field )
	{
	        flux = q*idt;	
           	fprintf( screen, "Old break down at the wall! (electricf= %e) (charge = %e)\n", electricf+electricf_ext, q  );
           	fprintf( screen, "Work = %e electric_local = %e electric = %e \n", dwork, perm*electricf, perm*electricf_ext );
	} 
   }

   dirFlux[0] = flux*delx;
   dirFlux[1] = flux*dely;
   dirFlux[2] = flux*delz;
   
   //Add half of the flux (located at the contact) to each particle in contact
   fix_efieldgran->efieldFlux[iPart] -= flux;
   fix_efieldgran->directionalEfieldFlux[iPart][0] -= 0.50 * dirFlux[0];
   fix_efieldgran->directionalEfieldFlux[iPart][1] -= 0.50 * dirFlux[1];
   fix_efieldgran->directionalEfieldFlux[iPart][2] -= 0.50 * dirFlux[2];
        
}


