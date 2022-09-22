/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Evan Smuts (U Cape Town, surface velocity rotation)
------------------------------------------------------------------------- */

#include "fix_mesh_surface_charge.h"
#include "fix_mesh_surface_charge_conduction.h"
#include <stdio.h>
#include "string.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "comm.h"
#include "math_extra.h"
#include "fix_property_global.h"
#include "fix_gravity.h"
#include "mpi.h"

#include "fix_charge_gran.h"
#include "math.h"

#include "efield_model.h"
#include "efield_model_normal.h"
#include "efield_model_screened.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMeshSurfaceChargeConduction::FixMeshSurfaceChargeConduction(LAMMPS *lmp, int narg, char **arg)
: FixMeshSurfaceCharge(lmp, narg, arg),
  fix_acc( NULL ),
  transfer_acceleration( 1.0 ),
  fix_work_model( NULL ),
  fix_work_a( NULL ),
  fix_work_b( NULL ),
  fix_work_wall( NULL ),
  fix_ef_coupling( NULL ),
  fix_resistivity_wall( NULL )
{
  this->chargeTransferActive = true;
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceChargeConduction::~FixMeshSurfaceChargeConduction()
{}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceChargeConduction::post_create()
{
    FixMeshSurfaceCharge::post_create();
    mesh()->prop().addElementProperty<ScalarContainer<double> >("normalElectricField","comm_exchange_borders","frame_invariant","restart_no"); 
}

void FixMeshSurfaceChargeConduction::init()
{
    FixMeshSurfaceCharge::init();
    int max_type = atom->get_properties()->max_type();
    
    fix_acc = static_cast<FixPropertyGlobal*>(modify->find_fix_property("acceleration","property/global","scalar",0,0,style,false));
    
    if( fix_acc ) fprintf( screen, "Acceleration factor: %e \n", fix_acc->compute_scalar() );
    
    fix_work_wall = static_cast<FixPropertyGlobal*>(modify->find_fix_property("work_wall","property/global","scalar",0,0,style));
    
    fix_work_model = static_cast<FixPropertyGlobal*>(modify->find_fix_property("workModel","property/global","peratomtype",max_type,0,style));
    fix_work_a = static_cast<FixPropertyGlobal*>(modify->find_fix_property("work_a","property/particle","peratomtype",max_type,0,style));
    fix_work_b = static_cast<FixPropertyGlobal*>(modify->find_fix_property("work_b","property/particle","peratomtype",max_type,0,style));
    
    normalElectricField_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("normalElectricField");
    normalElectricField_->setAll( 0.0 );
    
    if ( fix_tribo->get_relaxation_model_flag()  )
	fix_resistivity_wall = static_cast<FixPropertyGlobal*>(modify->find_fix_property("resistivity_wall","property/global","scalar",0,0,style));
   
    //check if there is CFD coupling for electric field
    fix_ef_coupling = static_cast<FixPropertyParticle*>(modify->find_fix_property("electricfield","property/particle","vector",3,0,style,false));
    
}


int FixMeshSurfaceChargeConduction::setmask()
{
    int mask = FixMeshSurfaceCharge::setmask();
    mask |= POST_FORCE;
    return mask;
}


// -- charge transfer --
void FixMeshSurfaceChargeConduction::compute_scalar_transport( LIGGGHTS::ContactModels::CollisionData& cdata, double * v_wall, int iMesh, int iTri )
{
    
    // -- force calculations --
    FixMeshSurfaceCharge::compute_scalar_transport( cdata, v_wall, iMesh, iTri );
    
    int ipart = cdata.i;
    if( !(atom->mask[ipart] & groupbit) ) return;
    
    int *type = atom->type;

    // -- do not calculate charge transfer for ghosts --
    if( ipart >= atom->nlocal ) return;
    
    const double chargePart = (fix_charge->vector_atom)[ipart];
    
    
    const double faceArea = calcArea( iTri );
    double idt = 1.0/update->dt; 
    
    const double permittivity = fix_tribo->permittivity;
    const double delta_charge = fix_tribo->delta_charge;
    const double bd_field     = fix_tribo->bd_field;    
    
    const double chargeFace = surfaceChargeDensity( iTri ) * faceArea;
    
    // -- effective radius of sphere that matches the face area -- (used for charge calculations)
    const double faceReff = sqrt( faceArea / (4 * M_PI ) ); 
    
    double** prev_loc = this->fix_prev_loc->array_atom;
    
    //fprintf( screen, "compute_scalar_transport (0) \n" ); //fixme
    /*
    if( !prev_loc ) fprintf( screen, "No previous location set!\n" );
    if( !fix_work_model ) fprintf( screen, "No fix_work_model set!\n" );
    if( !fix_work_a ) fprintf( screen, "No fix_work_a set!\n" );
    if( !fix_work_b ) fprintf( screen, "No fix_work_b set!\n" );
    if( !fix_tribo ) fprintf( screen, "No fix_tribo set!\n" );
    if( !fix_tribo->efieldModel() ) fprintf( screen, "No efieldModel() set!\n" );
    
    fprintf( screen, "compute_scalar_transport (1) \n" ); //fixme
    */
        
    const double* x = atom->x[ipart];
    double xMesh[3] = {0,0,0};
    
    double dist = calcNearestPoint( x, xMesh, iTri );
    
    //triMesh_->center( iTri, xMesh );
    
    const double r = atom->radius[ipart];

    // -- no contact between the mesh element and the face --
    if( dist > r )
    {
        return; 
    }
    
    //fprintf( screen, "compute_scalar_transport (2) \n" ); //fixme
    
    if( fix_acc ) transfer_acceleration = fix_acc->compute_scalar();
    const double chiStiff = 1.0/lmp->force->chiStiffnessScaling();    
    
    // -- calculate old and present overlap areas --
        
    double normalVector[3] = {0,0,0};
    double dist_old = calcNearestPoint( prev_loc[ipart], normalVector, iTri ); // -- use normal vector as dummy --
        
    for(int k = 0; k < 3; ++k )
    {
	normalVector[k] = x[k] - xMesh[k];
    }
    
    //fprintf( screen, "compute_scalar_transport (3) \n" ); //fixme

        
    //fprintf( screen, "dist = %e    %e \n", dist, dist_old );	
	
    for( int k = 0; k < 3; ++k )
       normalVector[k] /= dist;
    
    
    double cA = M_PI * r * relu(r - dist);    
    double cA_old = M_PI * r * relu(r-dist_old);
    
    //fprintf( screen, "compute_scalar_transport (4) \n" ); //fixme
    
    //fprintf( screen, "areas = %e    %e \n", cA, cA_old );
     
    if( cA_old < 0 ) //no collision on the last step
	cA_old = 0;   
     
    // -- change in the surface area --	
    double dA = cA - cA_old;	
    
    //fprintf( screen, "compute_scalar_transport (5) \n" ); //fixme
    
    double flux;
    /*
    //Neglect charge transfer for separation	
    if( dA < 0 ) 
    {
       return;
    }
    */
    

    
    //polarization contribution
    
    double polField = 0;
    // -- chashes before this --
    
    //fprintf( screen, "compute_scalar_transport (6) \n" ); //fixme
    
    
    if( fix_p )
    {
       double polarizationEfield[3] = {0,0,0};
       fix_tribo->efieldModel()->computeEfieldPolarization( normalVector, r, fix_p->array_atom[ipart], polarizationEfield );
       
       for( int k = 0; k < 3; ++k )
       {
           polField += polarizationEfield[k] * normalVector[k];
	   
	   //std::cout<<"polarizationEfield[k] is "<<polarizationEfield[k]<<std::endl;
	   //std::cout<<"normalVector[k] is "<<normalVector[k]<<std::endl;
       }
	   
    }
    
    else
    {
    	//std::cout<<"fix_p is null"<<std::endl;
    }
    
    //std::cout<<"polField is "<<polField<<std::endl;
 
     // -- normal "slowly varying" electric field at the face -- (OpenFoam normal vectors are outward, hence the minus)
    //const double ef = -normalElectricField( iTri );
    double ef = 0;
 
     
    if( fix_ef_coupling )
    {
       double* pef = fix_ef_coupling->array_atom[ipart];
       
       for( int k = 0; k < 3; ++k )
           ef += pef[k] * normalVector[k];       
        
    }
    
    //fprintf( screen, "compute_scalar_transport (6.5)  %d  \n", type[ipart]-1 );
    
   
        
    int charging_model_type = fix_work_model->compute_vector( type[ipart]-1 );
    
    //fprintf( screen, " charging_model_type = %d \n", charging_model_type );
    
    double workfPart = FixEfieldGran::work( 
    						r, 
    						fix_work_a->compute_vector( type[ipart]-1 ), 
    						fix_work_b->compute_vector( type[ipart]-1 ), 
						charging_model_type 
					  );
    
    //fprintf( screen, "compute_scalar_transport (7) \n" ); //fixme
    
    //double workWall = EV2JOULE * fix_work_wall->compute_scalar();
    
    double workWall = FixEfieldGran::work( 
    						faceReff, 
    						fix_work_wall->compute_scalar(), 
    						0, 
						WALL 
					  );
    
    const double dwork = permittivity / (delta_charge * ELECTRON_CHARGE ) * (workfPart - workWall);
    
    // -- electric field contribution coming from the face and particles --
    double electricf = fix_tribo->efieldModel()->implElectricField( chargePart, chargeFace, r, faceReff );

    if (dA > 0)
    {    
    	flux = chiStiff*transfer_acceleration * idt * dA * 					 //area difference
			     			( dwork - 					 //scaled workfunction difference
				          	  permittivity * ( electricf + ef + polField ) );//ambient electric field contribution
	
	//fprintf( screen, "before relaxation: %e", flux);
	
	if ( fix_tribo->get_relaxation_model_flag()  ) 
	{
		double part_resistivity = fix_tribo->resistivity[type[ipart]-1];
		double wall_resistivity = fix_resistivity_wall -> compute_scalar();
                double effective_resistivity = part_resistivity > wall_resistivity? part_resistivity: wall_resistivity;
                flux -= chiStiff*transfer_acceleration*idt*(electricf + ef)/effective_resistivity*cA;
	//fprintf( screen, "after relaxation: %e\n", flux);
	}
    }
    else if (dA < 0)
    {
	flux = 0.0;
	if ( fix_tribo->get_relaxation_model_flag()  )
        {
                double part_resistivity = fix_tribo->resistivity[type[ipart]-1];
                double wall_resistivity = fix_resistivity_wall -> compute_scalar();
                double  effective_resistivity = part_resistivity > wall_resistivity? part_resistivity: wall_resistivity;
                flux -= chiStiff*transfer_acceleration*idt*(electricf + ef)/effective_resistivity*cA;
        }	
    }
    else
    {
	flux = 0.0;
    }

    double impl_coff = 1.0 + chiStiff * transfer_acceleration * dA* fix_tribo->efieldModel()->implCoeff( r, faceReff );
   

    //if (  fix_tribo->get_relaxation_model_flag()  ) 
    
    // fprintf( screen, "dwork is %e, electricf is %e, ef is %e \n", dwork, electricf, ef ); 
    
    //fprintf( screen, "dA = %e  coff = %e  \n", dA, fix_tribo->efieldModel()->implCoeff( r, faceReff ) );
    
    flux /= impl_coff;
    
    //fprintf( screen, "impl_coeff = %e   flux = %e  \n", impl_coff, flux/idt ); //fixme
    
    
    if( bd_model_flag() )
    {
	//fprintf( screen, "Interior is using new breakdown model...");

	if( abs( electricf + ef + polField ) > bd_field)
	{
	    double surface_area_part = 4.0*M_PI*r*r;
	    double surface_area_face = faceArea;

	    if ( (electricf + ef + polField) >= 0 )
	    {
	    	    double q_12 = (bd_field - (electricf + ef + polField))/((surface_area_part + surface_area_face)/(surface_area_part*surface_area_face))*permittivity;
		    flux = q_12*idt;
	    }

	    else
	    {
	    	    double q_12 = -(bd_field + (electricf + ef + polField))/((surface_area_part + surface_area_face)/(surface_area_part*surface_area_face))*permittivity;
		    flux = q_12*idt;
	    }

	    fprintf( screen, "New break down in the wall! ( %e %e %e) (charge = %e, %e) \n", electricf, ef, polField ,chargePart, chargeFace );  	
	}		
    }
    
    
    double dirFlux[3] = {0,0,0};
    
    const double delx = x[0] - xMesh[0];
    const double dely = x[1] - xMesh[1];
    const double delz = x[2] - xMesh[2];
    
    //fprintf( screen, "compute_scalar_transport (8) \n" ); //fixme
    
    dirFlux[0] = flux*delx;
    dirFlux[1] = flux*dely;
    dirFlux[2] = flux*delz; 
    
    //fprintf( screen, "compute_scalar_transport (9) \n" ); //fixme
    
    // -- add half the flux in the face (given in terms of total charge) --
    this->surfaceChargeDensityFlux(iTri) += flux / ( faceArea * idt );
    
    //fprintf( screen, " iTri = %d iPart = %d chargePart = %e chargeFace = %e \n ", iTri, ipart, chargePart, chargeFace );
    
    //Add half of the flux (located at the contact) to particle in contact
    this->fix_efieldFlux->vector_atom[ipart] -= flux;
    this->fix_directionalEfieldFlux->array_atom[ipart][0] -= dirFlux[0];
    this->fix_directionalEfieldFlux->array_atom[ipart][1] -= dirFlux[1];
    this->fix_directionalEfieldFlux->array_atom[ipart][2] -= dirFlux[2];
    
}























