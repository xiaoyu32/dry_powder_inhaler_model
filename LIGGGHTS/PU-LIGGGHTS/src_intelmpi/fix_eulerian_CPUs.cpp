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


#include "cfd_datacoupling.h"
#include "fix_cfd_coupling.h"
#include "library_cfd_coupling.h"
#include "fix_property_particle.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "mpi.h"
#include "fix_eulerian_CPUs.h"
#include <iostream>
#include <fstream>
#include <vector>


//#include "polyhedron_search.h" /* deprecated algorithm (JK) */
#include "ray_search.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEulerianCPUs::FixEulerianCPUs(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fix_prevCPU( NULL )
{    

  debug_flag = false;
  ncpus_ = 1;
  n_cfd = 0;
  n_dem = 0;
  
  if (strcmp(style,"cpus") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix cpus command");

  int iarg = 3;
  bool hasargs = true;
  
  compress_flag = true;
  polyhedron_flag = true;
  
  //face merging tolerance
  tol = 1.0e-10;
  
  while(iarg < narg && hasargs)
  {
  
    hasargs = false;
    
    if (strcmp(arg[iarg],"debug") == 0) {
    
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'debug'");
      if(strcmp(arg[iarg+1],"yes") == 0)
      debug_flag = true;
      else if(strcmp(arg[iarg+1],"no") == 0)
      debug_flag = false;
      else error->all(FLERR,"Illegal debug option called");
      
      iarg +=2;
      hasargs = true; 
         
    }else if( strcmp(arg[iarg],"compress") == 0 ){ 
      
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'compress'");
      if(strcmp(arg[iarg+1],"yes") == 0)
      compress_flag = true;
      else if(strcmp(arg[iarg+1],"no") == 0)
      compress_flag = false;
      else error->all(FLERR,"Illegal debug option called");
      
      iarg +=2;
      hasargs = true;  
       
    }else if( strcmp(arg[iarg],"polyhedron") == 0 ){ 
      
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'polyhedron'");
      if(strcmp(arg[iarg+1],"yes") == 0)
      polyhedron_flag = true;
      else if(strcmp(arg[iarg+1],"no") == 0)
      polyhedron_flag = false;
      else error->all(FLERR,"Illegal debug option called");
      
      iarg +=2;
      hasargs = true;  
       
    }else if( strcmp(arg[iarg],"tolerance") == 0 ){ 
      
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'tolerance'");
      
      tol = atof( arg[iarg+1] );
      
      iarg +=2;
      hasargs = true;  
       
    }else{ 
      	error->fix_error(FLERR,this,"unknown keyword");
    }
  }
  
  int me,nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  
  localtagCFD_ = new int[nprocs]; 
  libraryTagBuffer_ = new int[nprocs];
  
  for( int i = 0; i < nprocs; ++i )
      localtagCFD_[i] = -1;
  
  std::ifstream inputPtr;
  inputPtr.open("../DEM/in.cpus");

  if (me == 0) {
    inputPtr >> ncpus_;
    fprintf(screen,"Opening ../DEM/in.cpus, nCpus = %d \n",ncpus_);
  }
  
  MPI_Bcast(&ncpus_, 1, MPI_DOUBLE, 0, world);
  coords_ = new double[6*ncpus_];
  
  //bounding box
  bounds__ = new double[6]; //FIXME: due to periodic boundaries, do something smarter
   
  if (me == 0) {
    for(int ii = 0; ii < 6*ncpus_; ++ii)
    {
        inputPtr >> coords_[ii];
    }
  }
  
  MPI_Bcast(coords_, 6*ncpus_, MPI_DOUBLE, 0, world);
   
  part_assigned = (bool*)malloc(PART_ASSIGNED_SIZE*sizeof(bool) ); 
  spart_assigned =  PART_ASSIGNED_SIZE;
  
  for( int i = 0; i < 3; ++i ){
      bounds__[2*i] = coords_[2*i];
      bounds__[2*i+1] = coords_[2*i+1];
  }
  
  for( int ic = 1; ic < ncpus_; ++ic )
  {
      for( int i = 0; i < 3; ++i )
      {  
	  if( coords_[6*ic+2*i] < bounds__[2*i] ) bounds__[2*i] = coords_[6*ic+2*i];
	  if( coords_[6*ic+2*i+1] > bounds__[2*i+1] ) bounds__[2*i+1] = coords_[6*ic+2*i+1]; 
      }
  }
  
  if( polyhedron_flag ) init_search( nprocs );
    
  if(debug_flag) fprintf(screen,"Debugging mpi_pu \n");
  
}

/* ---------------------------------------------------------------------- */

FixEulerianCPUs::~FixEulerianCPUs()
{

    if( part_assigned ) free( part_assigned );
    if( coords_ ) delete[] coords_;
    if( localtagCFD_ ) delete[] localtagCFD_;
    if( bounds__ ) delete[] bounds__;
    if( poly_search ) delete[] poly_search;
    
    if( libraryTagBuffer_ ) delete[] libraryTagBuffer_;
    
    part_assigned = NULL;
    coords_ = NULL;
    localtagCFD_ = NULL;
    bounds__ = NULL;
    poly_search = NULL;
    
}

void FixEulerianCPUs::init_search( int nprocs )
{ 
    
    int me;
    double* v;
    double buffer[100];
    
    MPI_Comm_rank(world,&me);
    
    poly_search = new RaySearch[nprocs]; // new PolyhedronSearch[nprocs]; /* deprecated algorithm (JK) */
        
    if( me == 0 )
    {
        
	//open file for reading face information
	ifstream file;
	fprintf(screen,"Opening in.cpus.faces! \n"); 
	file.open( "../DEM/in.cpus.faces" );
	
	std::string line;
	int proc = -1;
	int index[nprocs];
	vector<double> parsed;
	vector<double*> face;
	
	while( getline( file, line ) )
	{
	    
	    if( line.size() >= 10 && strncmp( "Processor:", line.c_str(), 9 ) == 0 ){
	        proc = parse_int( line );
		
		if( proc >= nprocs ){
		   error->fix_error(FLERR,this,"Processor number incompatible!");
		}
		
		if( proc < 0 ){
		   fprintf( screen, line.c_str() );
		   error->fix_error(FLERR,this,"Wrong processor number!");
		}
		index[proc] = 0;
		
		//set memory fre flag to true
		poly_search[proc].free_memory = 1;
		
	    }else if( line.size() >= 7 && strncmp( "Normal:", line.c_str(), 7 ) == 0 )
	    {
	        
		if( proc < 0 )
		{
		    fprintf( screen, line.c_str() );
		    error->fix_error(FLERR,this,"Wrong Processor Number!!");
		}
		
		parsed = parse_vector( line, 7 );
		
		if( parsed.size() < 3 )
		{
		    error->fix_error(FLERR,this,"Normal vector is missing components!");
		}
		
		v = new double[3];
				
		for( int j = 0; j < 3; ++j )
		    v[j] = parsed.at(j);

		parsed.clear();
		++index[proc];
		
		getline( file, line );
		
		parsed = parse_vector( line, 0 ); 
	       	get_face_vector( face, parsed );

		//add the face information to PolySearch
		poly_search[proc].add_face( face, v );
		parsed.clear();
		
	    }
	    	   
	} 
		
	//end of file io
	file.close();
 	
	vector<double*> face_normals;
	vector< vector<double* > > faces;
	
	//send face information to remaining processors
	for( int id = 1; id < nprocs; ++id )
	{
	    
	    for( proc = 0; proc <nprocs; ++proc )
	    {
	        
		face_normals = poly_search[proc].get_face_normals();
		faces = poly_search[proc].get_faces();
		
		MPI_Send( &proc, 1, MPI_INT, id, 99, MPI_COMM_WORLD );
		MPI_Send( &index[proc], 1,MPI_INT, id, 99, MPI_COMM_WORLD );
		
		for( int i = 0; i < index[proc]; ++i ){
		   
		    v = face_normals.at(i);
		    int nsend( 3*faces.at(i).size() );
		    MPI_Send( v, 3, MPI_DOUBLE, id, 99, MPI_COMM_WORLD ); 
		    MPI_Send( &nsend, 1, MPI_INT, id, 99, MPI_COMM_WORLD );
		    
		    for( int j = 0; j < faces.at(i).size(); ++j )
		    {
		         for( int k = 0; k < 3; ++k )
		             buffer[3*j+k] = faces.at(i).at(j)[k]; 
		    }
		    
		    MPI_Send( buffer, 3*faces.at(i).size(), MPI_DOUBLE, id, 99, MPI_COMM_WORLD );
		   
		}
		
		
	    }
	    
	}	
		
		
    }else{
    
        //receive boundary face information from root
        int proc;
	int nindex;
	int nrecv;
	double * u;
	MPI_Status status;
	vector<double*> face;
	
	for( int iproc = 0; iproc < nprocs; ++iproc )
	{
	     
	    MPI_Recv( &proc, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status );
	    MPI_Recv( &nindex, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status ); 
	    
	    //set memory fre flag to true
	    poly_search[iproc].free_memory = 1; 
	     
	    for( int i = 0; i < nindex; ++i ){
	        
		v = new double[3];
		MPI_Recv( v, 3, MPI_DOUBLE, 0 , 99, MPI_COMM_WORLD, &status );
		MPI_Recv( &nrecv, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status );
		MPI_Recv( buffer, nrecv, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status );
		
		face.clear();
		for( int j = 0; j < nrecv/3; ++j )
		{
		    u = new double[3];
		    for( int k = 0; k < 3; ++k )
		        u[k] = buffer[3*j+k];
			
		     face.push_back( u );
		}
		
		poly_search[proc].add_face( face, v );
		
	    } 
	     
	}
	
    }
    
    if( compress_flag ) 
    {
       if( me == 0 ) fprintf( screen, "Compressing face information...\n" );

       //compress faces (TODO:change the compression order)
       for( int iproc = 0; iproc < nprocs; ++iproc )
       {
           if( me == 0 ){ poly_search[iproc].print_flag = true; fprintf( screen, "Working on CPU %d... \n", iproc ); }
           poly_search[iproc].compress( tol );

       }
       
       if( me == 0 ) fprintf( screen, "Face information compressed.\n" );
    }
    
    for( int iproc = 0; iproc < nprocs; ++iproc )
        poly_search[iproc].initSearchTrees();
	
    if( me == 0 ){    
    
       //Sanity checking
       double xx[nprocs][3];
       bool foo;	
	
       for( int i = 0; i < nprocs; ++i )
       {
	 for( int j = 0; j < 3; ++j ){
             xx[i][j] = ( coords_[6*i+2*j] + coords_[6*i+2*j+1] )/2.0;
	 }
	 
	 
	 foo = poly_search[i].is_inside( xx[i] );
	 //double myPoint[3] = {9.231005e-03, -1.141453e-02, 4.643331e-03};

	 //foo = poly_search[i].is_inside( myPoint );
	 
	 if( !foo ){
	     fprintf( screen, "CPU %d failed sanity check: in.cpus and in.cpus.faces are incompatible! \n", i );
	     //fprintf( screen, "Test point: %f %f %f! \n",myPoint[0], myPoint[1], myPoint[2] );
	     fprintf( screen, "Test point: %f %f %f! \n", xx[i][0], xx[i][1], xx[i][2] );   
	     //poly_search[i].printFaces();
	 }
	
	  
       }
          
    }
    
    //MPI_Barrier( MPI_COMM_WORLD );    
    //error->fix_error(FLERR,this,"Stop~!!!!!");	
	
    /*if( me == 0 )
    {
    double y[3];
	
	y[0] = 1.098349e-02; 
        y[1] = -5.647646e-04;
	y[2] = 7.640144e-01;
       
        for( int i = 0; i < nprocs; ++i )
	{
	   bool bar = poly_search[i].is_inside( y );
	   
           cout<< "Processor: "<< i<< " =: " << bar << endl;
        }
    }
    
    MPI_Barrier( MPI_COMM_WORLD );
    error->fix_error(FLERR,this,"Stop~!!!!!");*/
    
} 


void FixEulerianCPUs::get_face_vector( vector<double*>& face, vector<double> parsed )
{
    
    int n = parsed.size()/3;
    double* v;
    face.clear();    
    
    for( int i = 0; i < n; ++i )
    {
        
	v = new double[3];
	
	for( int j = 0; j < 3; ++j )
	   v[j] = parsed.at( 3*i+j );
	   
	face.push_back( v );
	
    }
        
}

vector<double> FixEulerianCPUs::parse_vector( std::string line, int i )
{
    
    
    int j = i;
    vector<double> parsed;
    bool flag = false;
    double dd;
        
    while( j < line.size() )
    {
        
	char c = line.c_str()[j];
        	
	if( is_digit( c ) ){
	
	    dd = parse_double( line, &j, flag );
	    
	    if( flag ){
	        parsed.push_back( dd );
	    }else{
	        error->fix_error(FLERR,this,"Parsing failed due to unknown number!");
	    }
	    
	}else{
	    ++j;
	}
    }
    
    return parsed;
    
}

int FixEulerianCPUs::parse_int( std::string line )
{
   
    int j = 0;
    std::string nbuffer;
    bool flag = false;
    
    
    while( j < line.size() )
    {
        char c = line.c_str()[j];
	
	if( is_digit( c ) && c != 'e' ){
	    
	   flag = true;
	   nbuffer.push_back( c );
	    
	}else if( flag ) break;
	
	++j;
       
    } 
   
    return atoi( nbuffer.c_str() );
   
}

double FixEulerianCPUs::parse_double( std::string line, int* i, bool& success )
{
    int j = i[0];
    int index = 0;
    std::string nbuffer;
    
    bool flag = false;
    char c;
    while( j < line.size() )
    {
        
	c = line.c_str()[j];
	
	if( is_digit( c ) ){
	   
 	    flag = true;
	    
	    nbuffer.push_back( c ); 
	      
	}else{
	    if( flag ){
	        ++index;
	        break;
	    }
	}
	
	if( c == ',' ) break;
		
	++j;
	++index;
    }
    
    if( flag )
        success = true;
    else
        success = false;
    
    i[0] += index;
        
    return atof( nbuffer.c_str() );
    
}

bool FixEulerianCPUs::is_digit( char c )
{

    if( c == '0' ){
	return true;    
    }else if( c == '1' ){
        return true;
    }else if( c == '2' ){
	return true;
    }else if( c == '3' ){
	return true;
    }else if( c == '4' ){
        return true;
    }else if( c == '5' ){
        return true;
    }else if( c == '6' ){
        return true;
    }else if( c == '7' ){
        return true;
    }else if( c == '8' ){
        return true;
    }else if( c == '9' ){
        return true;
    }else if( c == '.' ){
        return true;
    }else if( c == '-' ){
	return true;    
    }else if( c == '+' ){
	return true;    
    }else if( c == 'e' ){
	return true;
    }else if( c == 'E' ){
	return true;
    }else{
        return false;
    }
    
}


void FixEulerianCPUs::post_create()
{

    fix_prevCPU = static_cast<FixPropertyParticle*>(modify->find_fix_property("fix_prevCPU","property/particle","scalar",1,0,this->style,false));
    if(!fix_prevCPU)
    {
      char* fixarg[9];
      fixarg[0]= (char *)"fix_prevCPU";
      fixarg[1]= (char *)"all";
      fixarg[2]= (char *)"property/particle";
      fixarg[3]= (char *)"fix_prevCPU";
      fixarg[4]= (char *)"scalar";
      fixarg[5]= (char *)"no";
      fixarg[6]= (char *)"yes";
      fixarg[7]= (char *)"no";
      fixarg[8]= (char *)"0.";
      fix_prevCPU = modify->add_fix_property_atom(9,fixarg,style);
    } 
  
}

void FixEulerianCPUs::init()
{
      fix_prevCPU = static_cast<FixPropertyParticle*>(modify->find_fix_property("fix_prevCPU","property/particle","scalar",0,0,style));
}

/* ---------------------------------------------------------------------- */

int FixEulerianCPUs::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

inline int FixEulerianCPUs::prevCPU( int i ) const
{
    int ic = fix_prevCPU->vector_atom[i];

    if( ic < 0 ) ic = 0;
    if( ic >= ncpus_ ) ic = ncpus_-1;
    return ic; 
}

void FixEulerianCPUs::end_of_step()
{

    int me;
    MPI_Comm_rank(world,&me);
    
    int couple_nevery_ = static_cast<FixCfdCoupling*>( locate_coupling_fix(lmp) )->couple_nevery_;
    
    //if( me == 0 ) fprintf( screen, "FixEulerianCPUs::end_of_step() (1)  \n" );
    
    
    FixCfdCoupling* fix_coupling = static_cast<FixCfdCoupling*>( locate_coupling_fix(lmp) );
    CfdDatacoupling* cfd_data_coupling = static_cast<CfdDatacoupling*>( fix_coupling->get_dc() );
    cfd_data_coupling->markUpdate( false );
    
    
    if( couple_nevery_ == 0 ){
        return;
    }
    
    
    int ts = update->ntimestep;
    if( ts % couple_nevery_ != 0 ) return;

    //if( me == 0 ) fprintf( screen, "LIGGGHTS time step = %d     couple_nevery_ = %d \n", ts, couple_nevery_ );
    //MPI_Barrier( MPI_COMM_WORLD );
    
    this->updateCpus();
 
    //if( me == 0 ) fprintf( screen, "FixEulerianCPUs::end_of_step() (2)  time = %d    step = %d  \n", ts, couple_nevery_ );

    
} 

void FixEulerianCPUs::updateCpus()
{

    int proc_id;
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
    if( proc_id == 0 ) fprintf( screen, "Updated CFD assignations \n" );

    cpu_indexes.resize( atom->nlocal ); //FIXME: correct this
    
    //sanity checking
    for( int i = 0; i < cpu_indexes.size(); ++i )
    {
        cpu_indexes[i] = -1; 
    }
    
    FixEulerianCPUs::reallocPartAssign( atom->nlocal );
    
    
    for( int ic = 0; ic < ncpus_; ++ic )
        localtagCFD_[ic]=0;
    
    FixEulerianCPUs::bin_atoms_init();  
          
    for(int ic = 0; ic < ncpus_; ic++)
    {
	FixEulerianCPUs::bin_atoms(ic);
    }
    
    int failedCounter = 0;
        
    for( int i = 0; i < atom->nlocal; ++i )
       if( !part_assigned[i] )
       { 
	  /* fprintf( 
	            screen, 
	   	    "Warning: No matching CFD domain found for particle x = [%e, %e, %e]. \n", 
	   	    atom->x[i][0], atom->x[i][1], atom->x[i][2] 
		   ); */
	 ++failedCounter;	   
         //int ic = fix_prevCPU->vector_atom[i];
	 //cpu_indexes[i] = ic;
	 //part_assigned[i] = true; 
	 //localtagCFD_[ic]++;	
	 
	 //fprintf( screen, "Warning: No matching CFD domain found for particle ID = %d  v = [ %e %e %e ]  \n", i, atom->v[i][0], atom->v[i][1], atom->v[i][2] );
	 		
       }
    
    if( failedCounter > 0 ) fprintf( screen, "Warning: %d/%d particles have left the CFD domain!\n", failedCounter, atom->nlocal );
    
    //if( proc_id == 0 ) fprintf( screen, "BEFORE liggghts_get_localtag \n" );   
       
    n_cfd = liggghts_get_localtag(lmp, false);
    n_dem = atom->nlocal;
    
    //if( proc_id == 0 ) fprintf( screen, "AFTER liggghts_get_localtag \n" );
    
    //sanity check
    //if( ntotal != n_dem ) printf( "Error! (2) \n" ); //error->fix_error(FLERR,this,"Number of Particles not Conserved! \n");

}


/* ----------------------------------------------------------------------*/

// -- try assigning particles to their previous CPUs --
void FixEulerianCPUs::bin_atoms_init()
{

    int i,ic,ibin;
    double **x = atom->x;
    int *mask = atom->mask;
    int nall = atom->nlocal;
    int *tag = atom->tag;
    int myrank;

    MPI_Comm_rank(world,&myrank);
 
    for (i=0; i<nall; i++)
    {
        ic = prevCPU(i);
			
        ibin = coord2bin(x[i], ic );
	
	if( ibin == 3 && !part_assigned[i] )
	{
	   cpu_indexes[i] = ic;
	   part_assigned[i] = true; 
	   localtagCFD_[ic]++;	
	}	
	
    }
    
}

void FixEulerianCPUs::bin_atoms(int ic)
{
  int i,ibin;
  double **x = atom->x;
  int *mask = atom->mask;
  int nall = atom->nlocal;
  int *tag = atom->tag;
  int myrank;
  
  MPI_Comm_rank(world,&myrank);
  
  for (i=0; i<nall; i++)
  {
      //if(! (mask[i] & groupbit)) continue;
      
      if(  prevCPU(i) == ic ) continue; // -- already checked --
      
      if( !part_assigned[i] ) ibin = coord2bin(x[i],ic);
                  
      if( ibin == 3 && !part_assigned[i] )
      {
         fix_prevCPU->vector_atom[i] = double(ic);
	 cpu_indexes[i] = ic;
	 part_assigned[i] = true; 
	 localtagCFD_[ic]++;	
      }

  }
  
} 

/* ----------------------------------------------------------------------
   map coord to grid
------------------------------------------------------------------------- */

int FixEulerianCPUs::coord2bin(const double *x, int ic) const
{
  int i,iCell;
  double iCell_lo,iCell_hi;
  double xx[3];
  
  
  iCell = 0;

  for(i=0;i<3;i++)
  {
      xx[i] = x[i];
      
      //FIXME: caused by periodic boundaries, do something different
      //periodic boundaries (handles the case when xx is outside the bounding box)
      if( xx[i] < bounds__[2*i] )   xx[i] = bounds__[2*i+1]-bounds__[2*i] + xx[i]; 
      if( xx[i] > bounds__[2*i+1] ) xx[i] = bounds__[2*i]-bounds__[2*i+1] + xx[i];  

      iCell_lo = (xx[i] - coords_[6*ic+2*i]);
      iCell_hi = (xx[i] - coords_[6*ic+2*i+1]);
      if( iCell_lo >= 0 && iCell_hi <= 0 )
	  iCell++;
  }
  
  //perform polyhedron search on the position
  if( iCell == 3 && polyhedron_flag ){
  
     if( !poly_search[ic].is_inside( xx ) ){	
         //fprintf( screen, "xx[%d]= %e %e %e \n", ic, xx[0], xx[1], xx[2] ); 
         iCell = 0;
     }
     
  }
  
  return iCell;
}


void FixEulerianCPUs::reallocPartAssign( int size )
{
    
    if( size > spart_assigned )
    {
        spart_assigned = size;
	part_assigned = (bool*)realloc( part_assigned, size*sizeof(bool) );
    }
    
    for( int i = 0; i < size; ++i )
    {
        part_assigned[i] = false;
    }
        
}



