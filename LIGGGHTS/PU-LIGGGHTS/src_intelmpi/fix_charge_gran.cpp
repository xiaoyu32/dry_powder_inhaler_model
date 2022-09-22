/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Efield
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


#include "fix_charge_gran.h"
#include "fix_scalar_transport_equation.h"
#include "group.h"
#include "stdlib.h"

#include "atom.h"
#include "compute_pair_gran_local.h"
#include "fix_property_particle.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "mech_param_gran.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair_gran.h"
#include "mpi_liggghts.h"

#include "efield_model.h"
#include "efield_model_normal.h"
#include "efield_model_screened.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEfieldGran::FixEfieldGran(LAMMPS *lmp, int narg, char **arg) : 
Fix(lmp, narg, arg),
new_bd_model_flag( false ),
relaxation_model_flag( false ),
efieldModel_( NULL )
{

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix efield/gran needs per particle radius and mass");

  //if(strcmp(arg[iarg++],"initial_charge"))
  //  error->fix_error(FLERR,this,"expecting keyword 'initial_charge'");
  //C0 = atof(arg[iarg++]);
    
  int iarg = 3;
  bool initial_charge_flag = false;
  //new_bd_model_flag = false;
  
  random_flag = false;
  initChargeOnParticles = false;
  
  flux_flag = 0;
  
  opposite_flag = false;
  
  delta_weights = NULL;
  delta_charges = NULL;
  delta_flag = false;
  number_of_delta = 0;
  
  constant_ef[0] = 0.0;
  constant_ef[1] = 0.0;
  constant_ef[2] = 0.0;
  
  multiply_charge = 1.0;
  multiply_flag = false;
  
  
  
  while( iarg < narg )
  {
      
      if( strcmp( arg[iarg], "external_field") == 0 )
      {

	  if( iarg+2 >= narg ) error->fix_error(FLERR,this,"Too few arguments for key word external_field!");
	  ++iarg;

	  constant_ef[0] = atof( arg[iarg] );
	  ++iarg;
	  constant_ef[1] = atof( arg[iarg] );
	  ++iarg;
	  constant_ef[2] = atof( arg[iarg] );

      }else if(strcmp(arg[iarg],"model") == 0)
      {
          ++iarg;
          if( iarg >= narg ){
	      fprintf( screen, "%s \n", arg[iarg-1] ); 
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'model'");
	  }
	  
	  ADDEFIELDMODEL(none,EfieldModel,efieldModel_)
	  ADDEFIELDMODEL(normal,EfieldModelNormal,efieldModel_)
	  ADDEFIELDMODEL(screened,EfieldModelScreened,efieldModel_)
	  
	  if( !efieldModel_ ) error->fix_error(FLERR,this,"Unknown electric field model. Available models are: 'none', 'normal' and 'screened'!");
	  	  
      }   
      
      if( strcmp( arg[iarg], "multiply") == 0 )
      {

	  if( iarg >= narg ) error->fix_error(FLERR,this,"Too few arguments for key word multiply!");
	  ++iarg;
	  
	  multiply_flag = true;
          multiply_charge = atof( arg[iarg] );
	  
      }

      
      if(strcmp(arg[iarg],"initial_charge") == 0)
      {
          
	  ++iarg;
          
	  if( iarg >= narg ){
	      fprintf( screen, "%s \n", arg[iarg-1] ); 
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'initial charge'");
	  }
	  
	  C0 = atof(arg[iarg]);
	  initial_charge_flag = true;
	  	  
      }
      
      
      if(strcmp(arg[iarg],"bd_model") == 0)
      {
      	 ++iarg;
	 
	 if( iarg >= narg ){
	      fprintf( screen, "%s \n", arg[iarg-1] ); 
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'bd_model'");
	 }
	 
	 int bd_model = atoi(arg[iarg]);
	 
	 
	 if (bd_model == 1){
	 	new_bd_model_flag = true;
		
	 } 
	 else{
	 	new_bd_model_flag = false;
	 }

	       
      }
     

      if(strcmp(arg[iarg],"relaxation_model") == 0)
      {  
         ++iarg;
         
         if( iarg >= narg ){
              fprintf( screen, "%s \n", arg[iarg-1] ); 
              error->fix_error(FLERR,this,"not enough arguments for keyword 'relaxation_model'");
         }
         
         //int relaxation_model = atoi(arg[iarg]);
         relaxation_model_flag = true;               
      }
      
      if( strcmp(arg[iarg],"initial_distribution") == 0)
      {
      	  
	  ++iarg;
	  
	  if( iarg >= narg ){
	      fprintf( screen, "%s \n", arg[iarg-1] ); 
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'initial distribution'");
	  }
	  
	  if( (narg-iarg)%2 != 0 )
	  {
	      error->fix_error(FLERR,this,"number of arguments for 'initial distribution' has to be even number");
	  }
	  
	  int inarg = (narg - iarg)/2;
	  	  
	  //delta function parameters	
	  if( delta_weights )
	  {
	     free( delta_weights );
	  }
	    
	  delta_weights = (double*)malloc( sizeof( double ) * inarg );
	  
	  if( delta_charges )
	  {
	      free( delta_charges ); 
	  }
	  
	  delta_charges = (double*)malloc( sizeof( double ) * inarg );
	  
	  
	  //fprintf( screen, "SCREEN: %d \n", inarg );
	  
	  //allocate parameters
	  for( int ii = 0; ii < inarg; ++ii )
	  {
	      
	      delta_weights[ii] = atof( arg[iarg] );
	      ++iarg;
	      
	      delta_charges[ii] = atof( arg[iarg] );
	      ++iarg;
	      
	  }
	  
	  delta_flag = true;
	  number_of_delta = inarg;
	  
	  if( iarg >= narg ) break;
	  
      }
      
      if(strcmp(arg[iarg],"opposite_charge") == 0)
      {
          ++iarg;
          if( iarg >= narg ){
	      fprintf( screen, "%s \n", arg[iarg-1] ); 
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'initial charge'");
	  }
	  
	  C0 = atof(arg[iarg]);
	  initial_charge_flag = true;
	  opposite_flag = true;	  
      }
      
      if(strcmp(arg[iarg],"random_charge") == 0)
      {
          ++iarg;
          if( iarg >= narg ){
	      fprintf( screen, "%s \n", arg[iarg-1] ); 
	      error->fix_error(FLERR,this,"not enough arguments for keyword 'initial charge'");
	  }
	  
	  C0 = atof(arg[iarg]);
	  initial_charge_flag = true;
	  random_flag = true;	  
      }
      
      if(strcmp(arg[iarg],"save_flux") == 0)
      {
          //start book keeping for the charge transfer
	  flux_flag = 1;
      }
      
      ++iarg;
      
  }
  
  if( !efieldModel_ ) error->fix_error(FLERR,this,"No electric field model specified: specify with keyword 'model'");  
  
  if( !initial_charge_flag )
     C0 = 0.0;
  
  pair_gran = NULL;

  fix_charge = fix_efieldFlux = fix_efieldSource = NULL;
  fix_ste = NULL;
  fix_directionalEfieldFlux = NULL;
  peratom_flag = 1;      
  size_peratom_cols = 0; 
  peratom_freq = 1;

  fix_ef_coupling = NULL;
  
  fix_work_a = NULL;
  work_a = NULL;
  
  fix_work_b = NULL;
  work_b = NULL;
	
  fix_work_model = NULL;
  work_model = NULL;

  fix_resistivity = NULL;
  resistivity = NULL;
  
  fix_prev_loc = NULL;
  prev_loc = NULL;
  
  fix_qflux = NULL;
  
  fix_electricf_r = NULL;
  
  scalar_flag = 1; 
  global_freq = 1; 

  cpl = NULL;

  FHG_init_flag = false;
  
  ef_coupling_flag = 0;

}

FixEfieldGran::~FixEfieldGran(){

  if( work_a ) delete [] work_a;
  work_a = NULL;
  
  if( work_b ) delete [] work_b;
  work_b = NULL;
  
  if( work_model ) delete [] work_model;
  work_model = NULL; 
  
  if( efieldModel_ ) delete efieldModel_;
  
  if( delta_charges )
  {
     free( delta_charges );
     delta_charges = NULL;
  }

  if( delta_weights )
  {
     free( delta_weights );
     delta_weights = NULL;
  }
   
}

/* ---------------------------------------------------------------------- */


void FixEfieldGran::post_create()
{
  
  
  fix_prev_loc = static_cast<FixPropertyParticle*>(modify->find_fix_property("prev_loc","property/particle","vector",3,0,this->style,false));
  if(!fix_prev_loc)
  {
    char* fixarg[11];
    fixarg[0]= (char *)"prev_loc";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"prev_loc";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"no";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"no";
    fixarg[8]= (char *)"0.";
    fixarg[9]= (char *)"0.";
    fixarg[10]= (char *)"0.";
    fix_prev_loc = modify->add_fix_property_atom(11,fixarg,style);
  } 
  
  fix_electricf_r = static_cast<FixPropertyParticle*>(modify->find_fix_property("fix_electricf_r","property/particle","vector",3,0,this->style,false));
  if(!fix_electricf_r)
  {
    char* fixarg[11];
    fixarg[0]= (char *)"fix_electricf_r";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"fix_electricf_r";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"no";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"yes";
    fixarg[8]= (char *)"0.";
    fixarg[9]= (char *)"0.";
    fixarg[10]= (char *)"0.";
    fix_electricf_r = modify->add_fix_property_atom(11,fixarg,style);
  } 
   
  // register directional flux
  fix_directionalEfieldFlux = static_cast<FixPropertyParticle*>(modify->find_fix_property("directionalEfieldFlux","property/particle","vector",3,0,this->style,false));
  if(!fix_directionalEfieldFlux)
  {
    char* fixarg[11];
    fixarg[0]= (char *)"directionalEfieldFlux";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"directionalEfieldFlux";
    fixarg[4]= (char *)"vector";
    fixarg[5]= (char *)"no";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"no";
    fixarg[8]= (char *)"0.";
    fixarg[9]= (char *)"0.";
    fixarg[10]= (char *)"0.";
    fix_directionalEfieldFlux = modify->add_fix_property_atom(11,fixarg,style);
  }
  
  fix_ste = modify->find_fix_scalar_transport_equation("efieldtransfer");
  //fix_ste = static_cast<FixChargeTransportEquation*>(modify->find_fix_id("ste_efieldtransfer"));
  if(!fix_ste)
  {
    char **newarg = new char*[15];
    newarg[0] = (char *) "ste_efieldtransfer";
    newarg[1] = group->names[igroup];
    newarg[2] = (char *) "transportequation/scalar";
    newarg[3] = (char *) "equation_id";
    newarg[4] = (char *) "efieldtransfer";
    newarg[5] = (char *) "quantity";
    newarg[6] = (char *) "charge";
    newarg[7] = (char *) "default_value";
    newarg[8] = new char[30];
    sprintf(newarg[8],"%e",C0);
    newarg[9] = (char *) "flux_quantity";
    newarg[10] = (char *) "efieldFlux";
    newarg[11] = (char *) "source_quantity";
    newarg[12] = (char *) "efieldSource";
    newarg[13] = (char *) "capacity_quantity";
    //newarg[14] = (char *) "thermalCapacity";
    newarg[14] = "none";
    modify->add_fix(15,newarg);

    delete [] newarg[8];
    delete [] newarg;
    
  }
  
  if( flux_flag == 0 ) return;
  
  fix_qflux = static_cast<FixPropertyParticle*>(modify->find_fix_property("qflux","property/particle","scalar",1,0,this->style,false));
  
  if( !fix_qflux )
  {
    char* fixarg[9];
    fixarg[0]= (char *)"qflux";
    fixarg[1]= (char *)"all";
    fixarg[2]= (char *)"property/particle";
    fixarg[3]= (char *)"qflux";
    fixarg[4]= (char *)"scalar";
    fixarg[5]= (char *)"no";
    fixarg[6]= (char *)"yes";
    fixarg[7]= (char *)"no";
    fixarg[8]= (char *)"0.";
    fix_qflux = modify->add_fix_property_atom(9,fixarg,style);
  }
  

}

/* ---------------------------------------------------------------------- */

void FixEfieldGran::updatePtrs(){

  charge = fix_charge->vector_atom;
  
  efieldFlux = fix_efieldFlux->vector_atom;
  efieldSource = fix_efieldSource->vector_atom;
  directionalEfieldFlux = fix_directionalEfieldFlux->array_atom;
  
  prev_loc = fix_prev_loc->array_atom;
  
}

/* ---------------------------------------------------------------------- */

double FixEfieldGran::copyCharge(int part)
{
   charge = fix_charge->vector_atom;
   return charge[part];
}

/* ---------------------------------------------------------------------- */

double FixEfieldGran::work(double r, double a, double b, int type)
{
   if( type == PARTICLE )
      return EV2JOULE * ( a + b/log(r*1000000) );   
   else if( type == WALL )
      return EV2JOULE * a;
   else
      return 0;
}

void FixEfieldGran::init(){
  
  //int max_type = pair_gran->mpg->max_type();
  int max_type = atom->get_properties()->max_type();
  prev_loc_flag = 0;
  	
  if ( work_a ) 
     delete []work_a;

  work_a = new double[max_type];
  
  fix_work_a = static_cast<FixPropertyGlobal*>(modify->find_fix_property("work_a","property/particle","peratomtype",max_type,0,style));
		
  if ( work_b ) 
     delete []work_b;
  
  work_b = new double[max_type];
  fix_work_b = static_cast<FixPropertyGlobal*>(modify->find_fix_property("work_b","property/particle","peratomtype",max_type,0,style));
  
   
  if( work_model ) 
     delete [] work_model;
  
  work_model = new int[max_type];
  fix_work_model = static_cast<FixPropertyGlobal*>(modify->find_fix_property("workModel","property/global","peratomtype",max_type,0,style));
  
  //set the medium permittivity and cutoff distance for charge transfer
  permittivity = (static_cast<FixPropertyGlobal*>(modify->find_fix_property("permittivity","property/global","scalar",0,0,style)))->compute_scalar();
  delta_charge = (static_cast<FixPropertyGlobal*>(modify->find_fix_property("deltaCharge","property/global","scalar",0,0,style)))->compute_scalar();
  bd_field =
(static_cast<FixPropertyGlobal*>(modify->find_fix_property("bd_field","property/global","scalar",0,0,style)))->compute_scalar();


  //Insert workfunction parameters
  for(int i=0;i< max_type; ++i)
  {
      work_a[i] = fix_work_a->compute_vector(i);
      work_b[i] = fix_work_b->compute_vector(i);
      work_model[i] = (int)fix_work_model->compute_vector(i);
  }

  if (get_relaxation_model_flag())
  {
  if (resistivity)
     delete [] resistivity;
   
   resistivity = new double[max_type];
   fix_resistivity = static_cast<FixPropertyGlobal*>(modify->find_fix_property("resistivity","property/particle","peratomtype",max_type,0,style));
   
   for (int i = 0; i < max_type; ++i)
   {
	resistivity[i] = fix_resistivity->compute_vector(i);
   }
  }
 
  if (!atom->radius_flag || !atom->rmass_flag) error->all(FLERR,"Please use a granular atom style for fix efield/gran");

  // check if a fix of this style already exists FIXME: return pointer to the object?
  //if(modify->n_fixes_style(style) > 1)
  //  error->fix_error(FLERR,this,"cannot have more than one fix of this style");

  if(!force->pair_match("gran", 0)) error->all(FLERR,"Please use a granular pair style for fix efield/gran");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  dnum = pair_gran->dnum_pair();
  dnum_mine = pair_gran->fix_extra_dnum_index(this);
  history_flag = pair_gran->is_history();

  //fix_ste = static_cast<FixChargeTransportEquation*>(modify->find_fix_id("ste_efieldtransfer"));
  fix_ste = modify->find_fix_scalar_transport_equation("efieldtransfer");
  if(!fix_ste) error->all(FLERR,"Fix efield/gran needs a fix transportequation/charge to work with");
    
  fix_charge = static_cast<FixPropertyParticle*>(modify->find_fix_property("charge","property/particle","scalar",0,0,style));
  fix_efieldFlux = static_cast<FixPropertyParticle*>(modify->find_fix_property("efieldFlux","property/particle","scalar",0,0,style));
  fix_efieldSource = static_cast<FixPropertyParticle*>(modify->find_fix_property("efieldSource","property/particle","scalar",0,0,style));
  fix_directionalEfieldFlux = static_cast<FixPropertyParticle*>(modify->find_fix_property("directionalEfieldFlux","property/particle","vector",0,0,style));
  fix_electricf_r = static_cast<FixPropertyParticle*>(modify->find_fix_property("fix_electricf_r","property/particle","vector",0,0,style));
  
  double total_charge = 0.0;
  double rnum = 0.0;
  
  efieldModel_->connectToProperties( force->registry );
  
  fix_init_charge = static_cast<FixPropertyGlobal*>(modify->find_fix_property("initial_charge","property/particle","peratomtype",max_type,0,this->style,false));
 
    
  //set initial charge
  if( random_flag )
  {
  
     for( int i = 0; i < atom->nlocal; ++i )
     { 
         rnum = C0 * random_double();
	 total_charge += rnum;
	 fix_charge->vector_atom[i] = rnum;
     }
     
     total_charge /= atom->nlocal;
     
     for( int i = 0; i < atom->nlocal; ++i )
     { 
	 fix_charge->vector_atom[i] -= total_charge;
     }
     
  }
  
  if( opposite_flag )
  {
     fprintf( screen, "Setting opposite charge %d ... \n", atom->nlocal );
     for( int i = 0; i < atom->nlocal; ++i )
     { 
	 fix_charge->vector_atom[i] = pow( -1.0, double(i) ) * C0;
     }

     
  }
  
  //draw random samples from mixture delta distribution
  if( delta_flag  )
  {
     
     for( int i = 0; i < atom->nlocal; ++i )
     {
        double p = random_double();
	double ww = 0;
	
	for( int j = 0; j < number_of_delta; ++j )
	{
	    
	    ww += delta_weights[j];
	    
	    if( p <= ww )
	    {
		fix_charge->vector_atom[i] = delta_charges[j];
		break;
	    }
	    

	}
     
     }     
     /*if( delta_charges )
     {
        free( delta_charges );
        delta_charges = NULL;
     }
     
     if( delta_weights )
     {
        free( delta_weights );
        delta_weights = NULL;
     }*/
     
  }
  
  if( multiply_flag )
  {

      fprintf( screen, "Multiplying charge %d ... \n", atom->nlocal );

      for( int i = 0; i < atom->nlocal; ++i )
      {
	  fix_charge->vector_atom[i] *= multiply_charge;
      }
  }
  
  
  
  //sets initial charge by particle type (overrides other qualifiers)
  if( fix_init_charge )
  {
  
      initChargeOnParticles = true;
        
      for( int i = 0; i < max_type; ++i )
      {
         fprintf( screen, "Charge for type: %d is %e. (%d) \n", (i+1), fix_init_charge->compute_vector( i ), atom->nlocal );
      }
          
  }else{
      
      fprintf( screen, "No initial charge per type specified. \n" );
      
  }
  
  updatePtrs();
  
  //check if there is CFD coupling for electric field
  fix_ef_coupling = static_cast<FixPropertyParticle*>(modify->find_fix_property("electricfield","property/particle","vector",3,0,style,false));
  
  if( fix_ef_coupling ) ef_coupling_flag = 1;
  else 			ef_coupling_flag = 0;
  
  
  
  
  FHG_init_flag = true;
  
}

/* ---------------------------------------------------------------------- */
/*
int FixEfieldGran::n_history_extra()
{   
    return 1;
}

bool FixEfieldGran::history_args(char** args)
{
    //provide names and newtonflags for each history value
    //newtonflag = 0 means that the value is same
    args[0] = (char *) "touchcharge";
    args[1] = (char *) "0";
    return true;
}
*/

/*
 Save particles current locations
*/
void FixEfieldGran::end_of_step(){
   
  //becomes TRUE after the first save
  prev_loc_flag = 1;
   
  updatePtrs();
   
  int * mask = atom->mask;
  int nlocal = atom->nlocal;
   
  for (int i = 0; i < nlocal; i++)
  {
     if (mask[i] & groupbit)
     {
        prev_loc[i][0] = atom->x[i][0];
        prev_loc[i][1] = atom->x[i][1];
        prev_loc[i][2] = atom->x[i][2];
     }
  }

  //update ghosts
  fix_prev_loc->do_forward_comm();
  
  
}

int FixEfieldGran::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEfieldGran::initial_integrate(int vflag)
{
  
  updatePtrs();
     
  //reset efield flux
  //sources are not reset
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  int i;
  
  for ( i = 0; i < nlocal; i++)
  {

     //efieldFlux[i] = 0;
     directionalEfieldFlux[i][0] = 0.;
     directionalEfieldFlux[i][1] = 0.;
     directionalEfieldFlux[i][2] = 0.;
     if( fix_qflux ) fix_qflux->vector_atom[i] = 0.0;
     
  }
  
  for( i = 0; i < nlocal; ++i )
  {
     fix_electricf_r->array_atom[i][0] = 0;
     fix_electricf_r->array_atom[i][1] = 0;
     fix_electricf_r->array_atom[i][2] = 0;
  }
  
  
  fix_electricf_r->do_forward_comm();
  
  //initialize particles only after injection  
  if( initChargeOnParticles && atom->natoms > 0 )
  {
      
      fprintf( screen, "Initialized charge on particles. \n" );
      
      initChargeOnParticles = false;
  
      int *type = atom->type;
          
      for( int i = 0; i < atom->nlocal; ++i )
      {
      	 fix_charge->vector_atom[i] = fix_init_charge->compute_vector( type[i]-1 );
      }
      
  }
  //fix_efieldFlux->do_forward_comm();
  
  //update ghosts
  //fix_efieldFlux->do_forward_comm();
  //fix_directionalEfieldFlux->do_forward_comm();
  
}

/* ---------------------------------------------------------------------- */

double FixEfieldGran::compute_scalar()
{
    return fix_ste->compute_scalar();

}

/* ---------------------------------------------------------------------- */

void FixEfieldGran::cpl_evaluate(class ComputePairGranLocal * cpl){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement cpl_evaluate().\n", mystyle);
  error->all(FLERR, emsg);

}

/* ---------------------------------------------------------------------- */

void FixEfieldGran::register_compute_pair_local(class ComputePairGranLocal *ptr){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement register_compute_pair_local().\n", mystyle);
  error->all(FLERR, emsg);

}

/* ---------------------------------------------------------------------- */

void FixEfieldGran::unregister_compute_pair_local(class ComputePairGranLocal *ptr){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement unregister_compute_pair_local().\n", mystyle);
  error->all(FLERR, emsg);

}

double FixEfieldGran::random_double()
{
    double a;

    a = double( rand() );

    a = 2.0*( a - RAND_MAX/2.0 )/RAND_MAX;

    return a;
}

double FixEfieldGran::abs( double a )
{
    if( a > 0 ) return a;
    
    return -a;
}
