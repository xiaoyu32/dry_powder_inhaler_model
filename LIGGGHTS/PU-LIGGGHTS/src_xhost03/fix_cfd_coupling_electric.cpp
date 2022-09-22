
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_electric.h"
#include "fix_property_atom.h"

#include "atom.h"
#include "compute_pair_gran_local.h"
#include "fix_property_particle.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "mech_param_gran.h"
#include "neigh_list.h"
#include "pair_gran.h"

#define M_PI 3.1415926
#define SMALL 1E-16

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingElectric::FixCfdCouplingElectric(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg)
{
    
    int iarg = 2;
    
    electricfield_flag = 1;
    
    electricfield_force_flag = 1;
    electricfieldGradientFlag = false;
    
    constant_ef[0] = 0.0;
    constant_ef[1] = 0.0;
    constant_ef[2] = 0.0;
    
    nullFlag[0] = 0;
    nullFlag[1] = 0;
    nullFlag[2] = 0;
    
    if (strcmp(this->style,"couple/cfd/electric") != 0) 
    {
        error->fix_error(FLERR,this,"unknown keyword");
    }
    
    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
	
	if( strcmp( arg[iarg], "electricfield") == 0 ){
	     
	     if( iarg >= narg ) error->fix_error(FLERR,this,"Too few arguments for key word electricfield!");
	     ++iarg;
	     
	     if( strcmp( arg[iarg], "yes" ) == 0 ) 	electricfield_flag = 1;
	     else 					electricfield_flag = 0;	
	     
	     //hasargs = false;
	}
	
	if( strcmp( arg[iarg], "set_zero") == 0 ){
	   
	   if( iarg >= narg ) error->fix_error(FLERR,this,"Too few arguments for key word set_zero!");
	   
	   ++iarg;
	   
	   bool increaseIarg = true;
	   
	   for( int i = 0; i < 3; ++i )
	   {	   	
	   
		increaseIarg = false;
		
		if( iarg >= (narg-1) ) break;
		
		if( strcmp( arg[iarg], "x" ) == 0 )
		{
	   	   increaseIarg = true;
		   nullFlag[0] = 1;
		}else if( strcmp( arg[iarg], "y" ) == 0 )
		{  
		   increaseIarg = true;
		   nullFlag[1] = 1;
		}else if( strcmp( arg[iarg], "z" ) == 0 )
		{ 
		   increaseIarg = true;
		   nullFlag[2] = 1;
	        }
		
		
		if( increaseIarg ) 	++iarg;
		else 			break;
		
	   }
	   
	}
	
	if( strcmp( arg[iarg], "electricfield_force") == 0 ){
	     
	     if( iarg >= narg ) error->fix_error(FLERR,this,"Too few arguments for key word electricfield!");
	     ++iarg;
	     
	     if( strcmp( arg[iarg], "yes" ) == 0 ) 	electricfield_force_flag = 1;
	     else 					electricfield_force_flag = 0;	
	     
	     //hasargs = false;
	}
	
	if( strcmp( arg[iarg], "external_field") == 0 )
	{
	    
	    if( iarg+2 >= narg ) error->fix_error(FLERR,this,"Too few arguments for key word external_field!");
	    ++iarg;
	    
	    constant_ef[0] = atof( arg[iarg] );
	    ++iarg;
	    constant_ef[1] = atof( arg[iarg] );
	    ++iarg;
	    constant_ef[2] = atof( arg[iarg] );
	    
	}
	
	if( strcmp( arg[iarg], "correction") == 0 ){
	     
	     if( iarg >= narg ) error->fix_error(FLERR,this,"Too few arguments for key word correction!");
	     ++iarg;
	     
	     if( strcmp( arg[iarg], "yes" ) == 0 ) 	electicfield_correction = 1;
	     else 					electicfield_correction = 0;	
	     
	}
	
	if( strcmp( arg[iarg], "gradient") == 0 )
	{
	     if( iarg >= narg ) error->fix_error(FLERR,this,"Too few arguments for key word gradient (yes/no)!");
	     ++iarg;	
	     
	     if( strcmp( arg[iarg], "yes" ) == 0 ) 	electricfieldGradientFlag = true;
	     else 					electricfieldGradientFlag = false;		     
	     	     
	}
	
	++iarg;
	
    }
    
    fix_p = NULL;
    fix_coupling = NULL;
    fix_charge = NULL;
    fix_electric_field = NULL;
    pair_gran = NULL;
    gran_flag = 0;
    
    fix_electric_field_gradient = NULL;
    
    ef_total[0] = 0.0;
    ef_total[1] = 0.0;
    ef_total[2] = 0.0;	      
	      
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingElectric::~FixCfdCouplingElectric()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingElectric::post_create()
{
     
     if( gran_flag == 0 ) return;
     
     fix_electric_field_gradient = static_cast<FixPropertyParticle*>(modify->find_fix_property("electricfield_gradient","property/particle","vector",9,0,this->style,false));
     
     if(!fix_electric_field_gradient)
     {
         char* fixarg[17];
         fixarg[0]= (char *)"electricfield_gradient";
         fixarg[1]= (char *)"all";
         fixarg[2]= (char *)"property/particle";
         fixarg[3]= (char *)"electricfield_gradient";
         fixarg[4]= (char *)"vector";
         fixarg[5]= (char *)"no";
         fixarg[6]= (char *)"yes";
         fixarg[7]= (char *)"no";
         
	 fixarg[8]= (char *)"0."; 	//xx
         fixarg[9]= (char *)"0."; 	//xy
         fixarg[10]= (char *)"0."; 	//xz ...
	 fixarg[11]= (char *)"0."; 
         fixarg[12]= (char *)"0.";
         fixarg[13]= (char *)"0.";
	 fixarg[14]= (char *)"0.";
         fixarg[15]= (char *)"0.";
         fixarg[16]= (char *)"0.";
	 
         fix_electric_field_gradient = modify->add_fix_property_atom(17,fixarg,style);
     }     
     

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingElectric::pre_delete(bool unfixflag)
{
    if(unfixflag && electricfield_flag) modify->delete_fix("electricfield");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingElectric::setmask()
{
    int mask = 0;
    mask |= POST_FORCE | PRE_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingElectric::init()
{
     
     
     
     if(modify->n_fixes_style(style) != 1)
         error->fix_error(FLERR,this,"More than one fix of this style is not allowed");
     
     //check if granular collision model is on, if on employ correction to ambient electric field to
     //avoid double counting of charges
     if( force->pair_match("gran", 0) ){ 
         pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
         gran_flag = 1;
     }else{
         gran_flag = 0;
     }

      // find coupling fix
     fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
     if(!fix_coupling) error->fix_error(FLERR,this,"Fix couple/cfd/electric needs a fix of type couple/cfd");
     
     fix_charge = static_cast<FixPropertyParticle*>(modify->find_fix_property("charge","property/particle","scalar",0,0,style));
     if(!fix_charge) error->fix_error(FLERR,this,"Fix couple/cfd/electric needs a fix property atom charge");
     
     //add charge to push properties
     fix_coupling->add_push_property("charge","scalar-atom");
     
     fix_p = static_cast<FixPropertyParticle*>(modify->find_fix_property("polarization","property/particle","vector",3,0,this->style,false));
     
     //add electric field to pull properties
     //if( electricfield_flag ) 
     
     permittivity = (static_cast<FixPropertyGlobal*>(modify->find_fix_property("permittivity","property/global","scalar",0,0,style)))->compute_scalar();
     
     fix_electric_field =
     static_cast<FixPropertyParticle*>(modify->find_fix_property("electricfield","property/particle","vector",3,0,style,false));
    
     if(!fix_electric_field)
     {
        const char* fixarg[11];
        fixarg[0]="electricfield";
        fixarg[1]="all";
        fixarg[2]="property/particle";
        fixarg[3]="electricfield";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="yes";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_electric_field = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
     }
	
     fix_electric_field_gradient = static_cast<FixPropertyParticle*>(modify->find_fix_property("electricfield_gradient","property/particle","vector",9,0,style,false));
     
     if(!fix_electric_field_gradient)
     {
         char* fixarg[17];
         fixarg[0]= (char *)"electricfield_gradient";
         fixarg[1]= (char *)"all";
         fixarg[2]= (char *)"property/particle";
         fixarg[3]= (char *)"electricfield_gradient";
         fixarg[4]= (char *)"vector";
         fixarg[5]= (char *)"yes";
         fixarg[6]= (char *)"no";
         fixarg[7]= (char *)"no";
         
	 fixarg[8]= (char *)"0."; 	//xx
         fixarg[9]= (char *)"0."; 	//xy
         fixarg[10]= (char *)"0."; 	//xz ...
	 fixarg[11]= (char *)"0."; 
         fixarg[12]= (char *)"0.";
         fixarg[13]= (char *)"0.";
	 fixarg[14]= (char *)"0.";
         fixarg[15]= (char *)"0.";
         fixarg[16]= (char *)"0.";
	 
         fix_electric_field_gradient = modify->add_fix_property_atom(17,fixarg,style);
     }


     
    
     //couple electricfield
     //fprintf(screen, "Adding Pull Property Electric field...\n");
     
     if( fix_coupling )
     {
	fix_coupling->add_pull_property("electricfield","vector-atom");
	fix_coupling->add_pull_property("electricfield_gradient","vector-atom");
     }
     //printf("Pull Property Electric field added\n");
     
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingElectric::pre_force(int)
{    
     //FixCfdCouplingElectric::compute_ef_correction();   
     
     int *mask = atom->mask;
     int nlocal = atom->nlocal;
          
     for( int ii = 0; ii < 3; ++ii )
     {
         
	 if( nullFlag[ii] == 1 )
	 {
	    for (int i = 0; i < nlocal; ++i)
	    {
	        if ( !(mask[i] & groupbit) ) continue;
	           fix_electric_field->array_atom[ii][i] = 0.0;
	    }
	 }
	 
     }
     
     fix_electric_field->do_forward_comm();
     
}

void FixCfdCouplingElectric::post_force(int)
{	
          
     double **f = atom->f;
     double** tf = atom->torque;
     double tor[3], tor2[3];
     
     int *mask = atom->mask;
     int nlocal = atom->nlocal;
     double **ef = fix_electric_field->array_atom;
     double *q = fix_charge->vector_atom;
     
     vectorZeroize3D(ef_total);
             
     // add electric force to force vector
     
     if( electricfield_force_flag == 1 )
     {

	for (int i = 0; i < nlocal; i++)
	{
            if (mask[i] & groupbit)
            {


		//add electric force
		for( int j = 0; j < 3; ++j )
		{
	           //fprintf( screen, "e= %e %e %e \n", ef[i][0], ef[i][1], ef[i][2] );
		   f[i][j] += q[i]*constant_ef[j];
		   f[i][j] += q[i]*ef[i][j];
		   ef_total[j] += ef[i][j];
		}  

		if( !fix_p ) continue;

		//polarization (only torque effects currently implemented)
		cross_product( fix_p->array_atom[i], constant_ef, tor );
		cross_product( fix_p->array_atom[i], ef[i], tor2 );

		//fprintf( screen, "p= %e %e %e \n", fix_p->array_atom[i][0], fix_p->array_atom[i][1], fix_p->array_atom[i][2] );
		//fprintf( screen, "e= %e %e %e \n", constant_ef[0], constant_ef[1], constant_ef[2] );
		//fprintf( screen, "t= %e %e %e \n", tor[0], tor[1], tor[2] );

		for( int j = 0; j < 3; ++j )
		{
	            tf[i][j] += tor[j] + tor2[j];
		}
		
		//force due to macroscopic electric field gradient
		if( electricfieldGradientFlag )
		{
		    
		    double fp[3] = {0,0,0};
		    double* gradE = fix_electric_field_gradient->array_atom[i];
		    
		    for( int ii = 0; ii < 3; ++ii )
		    	for( int jj = 0; jj < 3; ++jj )
			{
			    fp[ii] += gradE[jj + 3*ii] * fix_p->array_atom[i][jj];
			}
		    
		    for( int ii = 0; ii < 3; ++ii )
		       	f[i][ii] += fp[ii];
		    
		}
		
		
		
            }
	}

     }
     
}

/* ----------------------------------------------------------------------
   return components of total force on fix group
------------------------------------------------------------------------- */

double FixCfdCouplingElectric::sign( double* a, int k ) const
{
	
   if( a[k] < 0 )
      return -1.0;
   else
      return 1.0;

}

double FixCfdCouplingElectric::sign( double a ) const
{
	
   if( a < 0 )
      return -1.0;
   else
      return 1.0;

}

double FixCfdCouplingElectric::abs( double a ) const
{
   
   if( a < 0 )
 	return -a;
   else
	return a;	

}

inline double FixCfdCouplingElectric::dot( double* a, double* b ) const
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double FixCfdCouplingElectric::compute_vector(int n)
{
  MPI_Sum_Vector(ef_total,3,world);
  return ef_total[n];
}

void FixCfdCouplingElectric::cross_product( const double* a, const double* b, double* c ) const
{

    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];

}

























