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

#include <cstdlib>

#include "fix_charge_gran.h"
#include "fix_charge_gran_conduction.h"

#include "math.h"
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
#include "comm.h"				

#include "efield_model.h"
#include "efield_model_normal.h"
#include "efield_model_screened.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEfieldGranCond::FixEfieldGranCond(class LAMMPS *lmp, int narg, char **arg) : FixEfieldGran::FixEfieldGran(lmp, narg, arg)
{
  
  int iarg = 5;
  area_correction_flag = 0;
  
  bool hasargs = true;
  
  transfer_acceleration = 1.0;
  
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
  
  /*
  if( bd_model_flag() )
  {
      fprintf( screen, "Contructor is still using new breakdown model...");
  }
  */
  
  polarization_flag = false;
  fix_p = NULL;
}

/* ---------------------------------------------------------------------- */

FixEfieldGranCond::~FixEfieldGranCond()
{
  if( work_a ) delete [] work_a;
  work_a = NULL;
  
  if( work_b ) delete [] work_b;
  work_b = NULL;
  
  if( work_model ) delete [] work_model;
  work_model = NULL; 
  
}

/* ---------------------------------------------------------------------- */

// post_create() of parent is fine

/* ---------------------------------------------------------------------- */

void FixEfieldGranCond::pre_delete(bool unfixflag)
{

  // tell cpl that this fix is deleted
  if(cpl && unfixflag) cpl->reference_deleted();

}

/* ---------------------------------------------------------------------- */

int FixEfieldGranCond::setmask()
{
  int mask = FixEfieldGran::setmask();
  //mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;  
  return mask;
  //return 0;
}

/* ---------------------------------------------------------------------- */

void FixEfieldGranCond::init()
{
  
  const double *Y, *nu, *Y_orig;
  double expo, Yeff_ij, Yeff_orig_ij, ratio;
  int max_type;

  Fix *ymo_fix;
  
  if (FHG_init_flag == false){
     FixEfieldGran::init();
  }

  // calculate efield transfer correction
  
  max_type = atom->get_properties()->max_type();
  
  ymo_fix = NULL;
  if(area_correction_flag)
  {
    ymo_fix = modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",0,0,style);

    if(force->pair_match("gran/hooke",0)) expo = 1.;
    else if(force->pair_match("gran/hertz",0)) expo = 2./3.;
    else error->fix_error(FLERR,this,"area correction could not identify the granular pair style you are using, supported are hooke and hertz types");

    Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
    nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
    Y_orig = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_values();

    // allocate a new array within youngsModulusOriginal
    static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->new_array(max_type,max_type);

    // feed deltan_ratio into this array
    for(int i = 1; i < max_type+1; i++)
    {
      for(int j = 1; j < max_type+1; j++)
      {
        Yeff_ij      = 1./((1.-pow(nu[i-1],2.))/Y[i-1]     +(1.-pow(nu[j-1],2.))/Y[j-1]);
        Yeff_orig_ij = 1./((1.-pow(nu[i-1],2.))/Y_orig[i-1]+(1.-pow(nu[j-1],2.))/Y_orig[j-1]);
        ratio = pow(Yeff_ij/Yeff_orig_ij,expo);
        
        static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->array_modify(i-1,j-1,ratio);
      }
    }

    // get reference to deltan_ratio
    deltan_ratio = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_array_modified();
  }

  fix_p = static_cast<FixPropertyParticle*>(modify->find_fix_property("polarization","property/particle","vector",3,0,this->style,false));
    
  
  if( fix_p ) polarization_flag = true;
  
  updatePtrs();
  
  //error checks on coarsegraining
  if(force->cg_active())
     error->cg(FLERR,this->style);
    
}

/* ---------------------------------------------------------------------- */

void FixEfieldGranCond::post_force(int vflag){

  //chargelate function for using touchflag or not
  if(history_flag == 0) post_force_eval<0>(vflag,0);
  if(history_flag == 1) post_force_eval<1>(vflag,0);

}

/* ---------------------------------------------------------------------- */

void FixEfieldGranCond::cpl_evaluate(ComputePairGranLocal *caller)
{
  if(caller != cpl) error->all(FLERR,"Illegal situation in FixEfieldGranCond::cpl_evaluate");
  if(history_flag == 0) post_force_eval<0>(0,1);
  if(history_flag == 1) post_force_eval<1>(0,1);
}

/* ---------------------------------------------------------------------- */

template <int HISTFLAG>
void FixEfieldGranCond::post_force_eval(int vflag,int cpl_flag)
{
  
  //double hc,contactArea;
  double cA, cA_old, r_old;
  double flux,dirFlux[3],flux_old;
  double dA,electricf,electricf_ext;
  //double deltav, deltaf;
  double workfi, workfj, dwork;
  
  /*
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r; //rinv,rsqinv,tcoi,tcoj;
  
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *touchfix,*allhist,**firsthist;
  */

  int newton_pair = force->newton_pair;
  int inum = pair_gran->list->inum;
  int * ilist = pair_gran->list->ilist;
  int * numneigh = pair_gran->list->numneigh;
  int ** firstneigh = pair_gran->list->firstneigh;
  int *type = atom->type;
  
  double *radius = atom->radius;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double idt = 1.0/pair_gran->get_dt();
  
  int charging_model_typei;
  int charging_model_typej;
  
  if (strcmp(force->pair_style,"hybrid")==0)
    error->warning(FLERR,"Fix efield/gran/particle implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0)
    error->warning(FLERR,"Fix efield/gran/particle implementation may not be valid for pair style hybrid/overlay");
  
  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
    
  FixPropertyGlobal* fix_acc = static_cast<FixPropertyGlobal*>(modify->find_fix_property("acceleration","property/global","scalar",0,0,style,false));
   
  if( fix_acc ) transfer_acceleration = fix_acc->compute_scalar();
  const double chiStiff = 1.0/lmp->force->chiStiffnessScaling(); 
  
  updatePtrs();
  
  if( polarization_flag )
  {
      fix_p->do_forward_comm();
  }
  
  //evaluate the external field coming from neighboring particles
  for( int ii = 0; ii < inum; ++ii )
  {
     int i = ilist[ii];
     double xtmp = x[i][0];
     double ytmp = x[i][1];
     double ztmp = x[i][2];
     double radi = radius[i];
     int jnum = numneigh[i];
     int * jlist = firstneigh[i];  
     
     for (int jj = 0; jj < jnum; jj++) {
	
	int j = jlist[jj];
        j &= NEIGHMASK;
	
	if (i < j) continue;
	
	double delx = xtmp - x[j][0];
	double dely = ytmp - x[j][1];
	double delz = ztmp - x[j][2];
	double radj = radius[j];
	double rsq = delx*delx + dely*dely + delz*delz;
	double radsum = radi + radj;
     	
	//vector pointing from particle i to particle j
	double pij[3];
	double normal[3];
	double rd = 0.0;

	for( int k = 0; k < 3; ++k )
	{
	    pij[k] = (x[j][k] - x[i][k]);  
	    normal[k] = (x[j][k] - x[i][k]);
	    rd += pij[k]*pij[k];
	}  

	rd = sqrt( rd ); 
	
	for( int k = 0; k < 3; ++k )
	{
	    normal[k] /= rd;
	}
	
	//polarization contribution
	if( polarization_flag )
	{

	   efieldModel_->computeEfieldPolarization( normal, rd, fix_p->array_atom[j], fix_electricf_r->array_atom[i] );
	   if( j < nlocal ) efieldModel_->computeEfieldPolarization( normal, rd, fix_p->array_atom[i], fix_electricf_r->array_atom[j] );   

	}
	
	//particles are not colliding (count electric field to external near field)
	if (rsq > radsum*radsum) 
	{
	    
	   //compute electric field 
	   efieldModel_->computeEfieldCharge( normal, rd, charge[j], fix_electricf_r->array_atom[i] );
	   if( j < nlocal ) efieldModel_->computeEfieldCharge( normal, rd, -charge[i], fix_electricf_r->array_atom[j] );

	}
	
     }
  }
  
  //fix_electricf_r->do_reverse_comm();
  
  //update ghost particles
  fix_electricf_r->do_forward_comm(); 
  fix_prev_loc->do_forward_comm();
  
  if( ef_coupling_flag != 0 )
  {
     fix_ef_coupling->do_forward_comm();
  }
  
  // loop over neighbors of my atoms
  for(int ii = 0; ii < inum; ii++) {
     
    //xtmp = x[i][0];
    //ytmp = x[i][1];
    //ztmp = x[i][2];

    //radi = radius[i];
    //jlist = firstneigh[i];
    //jnum = numneigh[i];

    int i = ilist[ii];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    double radi = radius[i];
    int jnum = numneigh[i];
    int * jlist = firstneigh[i];   
        
    for (int jj = 0; jj < jnum; jj++) {
    	  
      //j = jlist[jj];
      //j &= NEIGHMASK;
      int j = jlist[jj];
      j &= NEIGHMASK;
      
      if (i < j) continue;
      
      if ( !(mask[i] & groupbit) && !(mask[j] & groupbit) ) continue;
      
      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double radj = radius[j];
      double rsq = delx*delx + dely*dely + delz*delz;
      double r = sqrt(rsq);
      double radsum = radi + radj;
 	
      //if ( true ) {	
      if ( rsq < radsum*radsum ) {  //contact
        
	/*
        if(HISTFLAG)
        {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radj = radius[j];
          radsum = radi + radj;
          if(rsq >= radsum*radsum) continue;
        }
	*/
        //r = sqrt(rsq);

        if(area_correction_flag)
        {
          //delta_n = radsum - r;
          //delta_n *= deltan_ratio[type[i]-1][type[j]-1];
          //r = radsum - delta_n;
        }
	
	
	//compute previous 
	if( prev_loc_flag ){
	    
	   r_old = 0.0;
	    
	   for(int k = 0; k < 3; ++k )
	       r_old += (prev_loc[i][k]-prev_loc[j][k])*(prev_loc[i][k]-prev_loc[j][k]);
	    
	   r_old = sqrt(r_old);
	   cA_old = M_PI * (radi*radj)/(radsum) * (radsum-r_old);
	   
	   if( cA_old < 0 ) //no collision on the last step
	       cA_old = 0;
	        
	}else{
	   cA_old = 0.0;
	}
	
	cA = M_PI * (radi*radj)/(radsum) * (radsum-r);    

	dA = cA - cA_old;


	//Neglect charge transfer for separating particles	
	/*
	if( dA < 0 )
	{  
	   dA = 0;
	   continue;
	}
	*/
	
	
	//ambient electric field contribution (assumes that delta << r)
	electricf_ext = 0.0;
	
	//unit vector pointing from particle j to particle i
	double pij[3];
	
	//vector pointing to "collision point"
	for( int k = 0; k < 3; ++k )
	{
	    //pc[k] = (x[j][k] + x[i][k])/2;
	    pij[k] = (x[i][k] - x[j][k])/r;  
	}
		
	//imposed electric field contribution
	for( int k = 0; k < 3; ++k )
	{
	    electricf_ext += pij[k] * constant_ef[k];
	}
		
	//contribution from surrounding particles
	for( int k = 0; k < 3; ++k )
	{
	    double eli, elj;    
	    
	    eli = fix_electricf_r->array_atom[i][k];
	    elj = fix_electricf_r->array_atom[j][k];
	    
	    electricf_ext += pij[k] * ( radj/radsum * eli + radi/radsum * elj );
	    
	    //if( j < nlocal ) electricf_ext += pij[k] * ( radj/radsum * eli + radi/radsum * elj );	
	    //else	     electricf_ext = pij[k]*eli;
	    
	}
	
	//fprintf( screen, "(1) flux=%e ele=%e \n", flux, electricf_ext );
	
	//electricf_ext /= (4*M_PI*permittivity);
	
	//external and ambient electricfield contribution
	if( ef_coupling_flag != 0 )
	{
	    
	    //linear interpolation for intersection point electricfield
	    for( int k = 0; k <3; ++k )
	    {
	    	//electricf_ext += pij[k] * ( 	radj/radsum * (fix_ef_coupling->array_atom[i][k] )+
	    	//               			radi/radsum * (fix_ef_coupling->array_atom[j][k] ) );
	    
	    	electricf_ext += pij[k] * ( 	radj/radsum * (fix_ef_coupling->array_atom[i][k] ) +
	    	               			radi/radsum * (fix_ef_coupling->array_atom[j][k] ) );
	    }
	    
	}
		
	//FIX: might need scaling, ELECTRON_CHARGE is very small -> Its ok
	
	charging_model_typei = work_model[type[i]-1];
	charging_model_typej = work_model[type[j]-1];
	//flux from i to j
	
	workfi = FixEfieldGran::work( radius[i], work_a[type[i]-1], work_b[type[i]-1], charging_model_typei );
	workfj = FixEfieldGran::work( radius[j], work_a[type[j]-1], work_b[type[j]-1], charging_model_typej );
	dwork = permittivity / (delta_charge*ELECTRON_CHARGE ) *(workfi - workfj);
        
	//int it = 0;
	
	flux = 0.0;
	flux_old = 1.0;
			
		
	//contribution from colliding particles
	//electricf = -charge[i]/(radi*radi) + charge[j]/(radj*radj);
	//electricf /= (4*M_PI*permittivity);
		
	electricf = efieldModel_->implElectricField( charge[i], charge[j], radi, radj ); //JK: should this be here?
	
	if (dA > 0)
	{
		flux = chiStiff*transfer_acceleration * idt * dA * 					//area difference
			     		( dwork - 						//scaled workfunction difference
				          permittivity*( electricf + electricf_ext ) ); 	//ambient electric field contribution
		//fprintf( screen, "before relaxation flux=%e", flux);
		//fprintf( screen, "dA: %e", dA);
		if (get_relaxation_model_flag())
		{
			double effective_resistivity = resistivity[type[i]-1] > resistivity[type[j]-1]? resistivity[type[i]-1]:resistivity[type[j]-1];
			flux -= chiStiff*transfer_acceleration*idt*(electricf + electricf_ext)/effective_resistivity*cA; 
			//fprintf( screen, "resistivity: %e", effective_resistivity);
			//fprintf( screen, "cA: %e", cA);
			//fprintf( screen, "after relaxation flux=%e\n", flux);       
		}
	}

	else if (dA < 0)
	{
		flux = 0.0;
		
 		if (get_relaxation_model_flag())
                {
                        double effective_resistivity = resistivity[type[i]-1] > resistivity[type[j]-1]? resistivity[type[i]-1]:resistivity[type[j]-1];
                        flux -= chiStiff*transfer_acceleration*idt*(electricf + electricf_ext)/effective_resistivity*cA;
                }
	}
	
	else
	{
		flux = 0.0;
	}	
	
	//fprintf( screen, "dwork is %e, permittivity is %e, electricf is %e, electricf_ext is %e \n", dwork, permittivity, electricf, electricf_ext );  
	
	//double impl_coff = 1.0 + chiStiff*transfer_acceleration*dA*( 1.0/(4*M_PI*radi*radi) + 1.0/(4*M_PI*radj*radj) );
	double impl_coff = 1.0 + chiStiff*transfer_acceleration*dA* efieldModel_->implCoeff( radi, radj ); 
	
	
	
	flux /= impl_coff;
	
	//contribution from colliding particles (after accounting the flux)
	//electricf = -(charge[i]-flux/idt)/(radi*radi) + (charge[j]+flux/idt)/(radj*radj);
	//electricf /= (4*M_PI*permittivity);
        //electricf = efieldModel_->implElectricField( charge[i]-flux/idt, charge[j]+flux/idt, radi, radj );


		

	if( bd_model_flag() )
	{
		
		if( abs( electricf + electricf_ext ) > bd_field)
		{
			double surface_area_i = 4.0*M_PI*radi*radi;
	    		double surface_area_j = 4.0*M_PI*radj*radj;
	    
	    		if ( (electricf + electricf_ext) >= 0)
	    		{
	    			double q_12 = (bd_field - (electricf + electricf_ext))/((surface_area_i + surface_area_j)/(surface_area_i*surface_area_j))*permittivity;
				//flux = 1.1*q_12*idt;
				flux = q_12*idt;
	    		}
	    
	    		else
	    		{
	    			double q_12 = -(bd_field + (electricf + electricf_ext))/((surface_area_i + surface_area_j)/(surface_area_i*surface_area_j))*permittivity;
				//flux = 1.1*q_12*idt;
				flux = q_12*idt;
	    		}
	    		
			fprintf( screen, "New break down in the interior! ( %e %e ) (charge = %e, %e) \n", electricf, electricf_ext, charge[i], charge[j] );  	
	 	}		
	} 
	
	else
	{
		
		if( abs( electricf + electricf_ext ) > bd_field)
		{
			double qsat = ( charge[i] + charge[j] )/2.0;

	    		flux = (charge[i] - qsat)*idt;

	    		fprintf( screen, "Old break down in the interior! (%e %e) ( charge= %e, %e ) \n", electricf, electricf_ext, charge[i], charge[j] );
		}
	}
	

	
        dirFlux[0] = flux*delx;
        dirFlux[1] = flux*dely;
        dirFlux[2] = flux*delz;
        
	if( !cpl_flag )
        {
	
	  //if( j >= nlocal ) continue;
	
          //Add half of the flux (located at the contact) to each particle in contact
          efieldFlux[i] -= flux;
          directionalEfieldFlux[i][0] -= dirFlux[0];
          directionalEfieldFlux[i][1] -= dirFlux[1];
          directionalEfieldFlux[i][2] -= dirFlux[2];
	  
	  //if( j >= nlocal ) continue;
	  if (newton_pair || j < nlocal)
	  {      
             efieldFlux[j] += flux;
             directionalEfieldFlux[j][0] +=  dirFlux[0];
             directionalEfieldFlux[j][1] +=  dirFlux[1];
             directionalEfieldFlux[j][2] +=  dirFlux[2];
          }
        }

        // if(cpl_flag && cpl) cpl->add_charge(i,j,flux); ERROR why ??
      }
    }
  }
  
  //fix_efieldFlux->do_reverse_comm();
  
  if(newton_pair) fix_efieldFlux->do_reverse_comm();
  //if(newton_pair) fix_directionalEfieldFlux->do_reverse_comm();
  
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void FixEfieldGranCond::register_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != NULL)
      error->all(FLERR,"Fix efield/gran/particle allows only one compute of type pair/local");
   cpl = ptr;
}

void FixEfieldGranCond::unregister_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != ptr)
       error->all(FLERR,"Illegal situation in FixEfieldGranCond::unregister_compute_pair_local");
   cpl = NULL;
}


//void FixEfieldGranCond::initial_integrate(int vflag)
//{
  /*updatePtrs();
     
  //reset efield flux
  //sources are not reset
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  int i;
  
  for ( i = 0; i < nlocal; i++)
  {
     if (mask[i] & groupbit)
     {
        efieldFlux[i] = 0;
        directionalEfieldFlux[i][0] = 0.;
        directionalEfieldFlux[i][1] = 0.;
        directionalEfieldFlux[i][2] = 0.;
	if( fix_qflux ) fix_qflux->vector_atom[i] = 0.0;
     }
  }
 
  //update ghosts
  //fix_efieldFlux->do_forward_comm();
  fix_directionalEfieldFlux->do_forward_comm();*/
  
//}



