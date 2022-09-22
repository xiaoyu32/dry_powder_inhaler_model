#include "wallComm.h"

#include "mpi.h"

#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "library_cfd_coupling.h"
#include "wallElementFinder.h"
          
using namespace LAMMPS_NS;

WallComm::WallComm() :
MPI_Communicator(),
n_cfd( 0 ),
n_dem( 0 ),
coords_( NULL ),
bounds__( NULL ),
wallElementFinder_( NULL ),
active( false ),
fileRead( false ),
n_send( NULL ),
n_recv( NULL )
{
   MPI_Comm_rank( MPI_COMM_WORLD, &id );
   MPI_Comm_size( MPI_COMM_WORLD, &ntask );

    // -- read CPU bounds --
   std::ifstream inputPtr;
   inputPtr.open("../DEM/in.cpus");
   
   if( inputPtr.is_open() )
   {  
       active = true;
       fileRead = true;
   }else{
       // -- no corresponding cpu information found, coupling is not in use --
       inputPtr.close();
       return;
   }
   
   
   if( id > 0  ) 
   {
       inputPtr.close();
   }
   
   int ncpus_ = 1;
    
   if (id == 0) {
     inputPtr >> ncpus_;
   }

   MPI_Bcast(&ncpus_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   coords_ = new double[6*ncpus_];

   // -- CFD task bounding box --
   bounds__ = new double[6]; //FIXME: due to periodic boundaries, do something smarter
   
   
   if (id == 0) {
            
     for(int ii = 0; ii < 6*ncpus_; ++ii)
     {
         inputPtr >> coords_[ii];
     }
     
     inputPtr.close();
     
   }
   
   
   MPI_Bcast(coords_, 6*ncpus_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   
   
   
   MPI_Barrier( MPI_COMM_WORLD );
   
   
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

   localtagCFD_ = new int[ntask];
   n_send = new int[ ntask ]; 
   n_recv = new int[ ntask ];

   for( int i = 0; i < ntask; ++i )
   {
      localtagCFD_[i] = 0;
      n_send[i] = 0;
      n_recv[i] = 0;
   } 

   elementCPUMapping_ = new std::vector<int>[ ntask ];   
   elementDEMMapping_ = new std::vector<int>[ ntask ]; 
   
   // -- initialize wall element finders --
   wallElementFinder_ = new WallElementFinder*[ntask];
   
   for( int itask = 0; itask < ntask; ++itask )
   {
       
       double lowBound[3] = {0,0,0};
       lowBound[0] = coords_[6*itask];
       lowBound[1] = coords_[6*itask+2];
       lowBound[2] = coords_[6*itask+4];
       
       double highBound[3] = {0,0,0};
       highBound[0] = coords_[6*itask+1];
       highBound[1] = coords_[6*itask+3];
       highBound[2] = coords_[6*itask+5];       
       
       wallElementFinder_[itask] = new WallElementFinder( lowBound, highBound );
   }
   
   
   if( id == 0 )
   {

       //open file for reading face information
       std::ifstream file;
       file.open( "../DEM/in.cpus.faces" );

       std::string line;
       int proc = -1;

       std::vector<double*> face;
       std::vector<double> parsingBuffer;
	
       int failedFaces = 0;	
	
       while( getline( file, line ) )
       {

	   if( line.size() >= 10 && strncmp( "Processor:", line.c_str(), 9 ) == 0 ){
	       proc = parse_int( line );

	       if( proc >= ntask ){
	          printf( "Error: incorrect processor number! \n" );
	       }

	       if( proc < 0 || proc >= ntask ){
		  printf( "Error: incorrect processor number! \n" );
	       }
	       
	   }else if( line.size() >= 7 && strncmp( "Normal:", line.c_str(), 7 ) == 0 )
	   {

	       if( proc < 0 )
	       {
		   printf( "Error: incorrect processor number! \n" );
	       }

	        
	       // -- do not read the normal vector --	       
	       getline( file, line );

	       parsingBuffer = parse_vector( line, 0 ); 
	       get_face_vector( face, parsingBuffer );

	       //add the face information to PolySearch
	       
	       if( !wallElementFinder_[proc]->addFace( face ) )
	       {
		  ++failedFaces; 
	       }
	       	       
	       for( int ii = 0; ii < face.size(); ++ii )
	       {  
	           delete[] face.at(ii);
	       }
	       
	       face.clear();
	       parsingBuffer.clear();

	   }

       } 

       //end of file io
       file.close();
       
       if( failedFaces > 0 ) printf( "Warning: failed to add %d CFD faces! \n", failedFaces );
       
       //send face information to remaining processors
       for( int itask = 0; itask < ntask; ++itask )
       {
           int ptrLength = wallElementFinder_[itask]->byteSize();
	   	   
	   char* elementFinderPtr = wallElementFinder_[itask]->extractByteStream();
	   
	   // -- number of bytes in the char presentation --
	   MPI_Bcast( &ptrLength, 1, MPI_INT, 0, MPI_COMM_WORLD );
	   
	   // -- sent the encoded elementFinder --
	   MPI_Bcast( elementFinderPtr, ptrLength, MPI_BYTE, 0, MPI_COMM_WORLD );
	  
	   delete[] elementFinderPtr;
	   
       }

   }else{
   
   	
       for( int itask = 0; itask < ntask; ++itask )
       {
           
	   int ptrLength = 0;
	   MPI_Bcast( &ptrLength, 1, MPI_INT, 0, MPI_COMM_WORLD );
	   
	   char* elementFinderPtr = new char[ ptrLength ];
	   
	   // -- receive the encoded elementFinder --
	   MPI_Bcast( elementFinderPtr, ptrLength, MPI_BYTE, 0, MPI_COMM_WORLD );
	   wallElementFinder_[itask]->initFinder( elementFinderPtr, ptrLength );
	   
	   delete[] elementFinderPtr;
       }	

   }
   
   //for( int itask = 0; itask < ntask; ++itask )
   //{
   //    printf( "ID = %d :  CFD = %d    number of elements = %d \n", id, itask, wallElementFinder_[itask]->numberOfElements() );   
   //} 
   
   MPI_Barrier( MPI_COMM_WORLD );
   
}

WallComm::~WallComm()
{
    // -- nothing to delete --
    if( !fileRead ) return;

    if( localtagCFD_ ) delete[] localtagCFD_;
    if( elementCPUMapping_ ) delete[] elementCPUMapping_;
    if( elementDEMMapping_ ) delete[] elementDEMMapping_;

    if( n_send ) delete[] n_send;
    if( n_recv ) delete[] n_recv;
    
    if( coords_ ) delete[] coords_;
    if( bounds__ ) delete[] bounds__;
    
    if( wallElementFinder_ )
    {
       for( int i = 0; i < ntask; ++i )
          delete wallElementFinder_[i];
       delete[] wallElementFinder_;
    }
    
}

// -- reconstruct MPI topology for further communication -- 
void WallComm::init( double** arr_, int nlocal )
{
    if( !fileRead ) return;    
    if( !isActive() ) return;
           
    // -- number of elements in the DEM --
    n_dem = nlocal;
    n_cfd = 0;
    
    cpu_indexes.clear();
    
    for( int ic = 0; ic < ntask; ++ic )
       localtagCFD_[ic] = 0;
    
    for( int i = 0; i < n_dem; ++i )
       cpu_indexes.push_back( -1 );
    
    // -- initialize number of elements sent and received from CFD part --   
    for( int itask = 0; itask < ntask; ++itask )
    {
        n_send[itask] = 0;
	n_recv[itask] = 0;
    }
    
    int multipleMatches = 0;
	
    for( int i = 0; i < n_dem; ++i )
    {
        
	const double* pos = arr_[i];
	
	int matchFound = 0;
	
	for( int ic = 0; ic < ntask; ++ic )
	{
	    
	    if( this->coord2bin( pos, ic ) == 3 ) // -- center lies in the bounding box --
	    {
	        		
		// 
		if( wallElementFinder_[ic]->isOnFaces( pos ) ) // -- check if center is on the collection of CFD boundary faces --
		{
		
		   if( matchFound == 0 )
		   {
		       cpu_indexes[i] = ic;
		      ++localtagCFD_[ic];
		   }

		   ++matchFound;		
		}
		//break; // -- FIXME: may need to send element to two tasks? --
	    }
	    
	    if( matchFound > 1 )
	    {
	        ++multipleMatches;
	    }
	    
	}
	
	/*if( cpu_indexes[i] < 0 )  // -- this should never happen --
	{
	    printf( "Warning: DEM wall element %d not sent to any processor! [ %e %e %e ] failure mode = %d \n", 
	    		i, arr_[i][0], arr_[i][1], arr_[i][2], cpu_indexes[i] );
	}*/
	
    }
    
    int failedElements = 0;
    
    for( int i = 0; i < n_dem; ++i )
    {
       if( cpu_indexes[i] < 0 ) ++failedElements;
    }
        
    if( multipleMatches > 0 ) printf( "Warning: %d DEM wall elements matched multiple CFD domains! (MPI task %d) \n", multipleMatches, id );
    if( failedElements > 0 )  printf( "Warning: Failed to find matching CFD domains for %d DEM wall elements! (MPI task %d) \n", failedElements, id );
	
    // -- compute the number of elements on the CFD side --
    int* buffer = new int[ntask];
    
    for( int i = 0; i < ntask; ++i )
        buffer[i] = 0;
    
    MPI_Allreduce( localtagCFD_, buffer, ntask, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    
    n_cfd = buffer[id];
    
    delete[] buffer;    
    
    // reconstruct element CPU mapping
    for( int ic = 0; ic < 0; ic < ntask )
       elementCPUMapping_[ic].clear();
    
    for( int i = 0; i < n_dem; ++i )
    {
       if( cpu_indexes[i] >= 0 ) elementCPUMapping_[ cpu_indexes[i] ].push_back( i ); // indexes that are sent to cpu
    }
        
    // reconstruct DEM mapping
    
    // clear former mapping
    for( int ic = 0; ic < 0; ic < ntask )
       elementDEMMapping_[ic].clear();    
    
    local_ids.clear();
    cpu_ids.clear();
        
    // fill local_ids and cpu_ids
    do_bcast_comm();
        
    // -- reconstruct DEM mapping --
    for( int i = 0; i < n_cfd; ++i )
       elementDEMMapping_[ cpu_ids[i] ].push_back( i );
 
  
    // -- mark the DEM side topology changed --
    bool& wallUpdated = isWallUpdated();
    wallUpdated = true;
    
    //printf( "ID = %d receiving = %d \n", id, n_cfd );
    
}

void WallComm::send_data( int dest_id )
{
    
    //init_send = true;    
    
    std::vector<int>& mapping = elementCPUMapping_[dest_id];
    
    n_send[dest_id] = mapping.size();//cpu_mapping->localtagCFD()[dest_id];
    
    int i;
    
    double * buff2;
    
    //send the number of particles that dest_id will be receiving
    MPI_Isend( &n_send[dest_id], 1, MPI_INT, dest_id, 99, MPI_COMM_WORLD, &destReqs ); 
     
    if( n_send[dest_id] == 0 )
    {
        //MPI_Wait( &send_reqs[0], &send_stats[0] );    
        //init_send = true;
	return;
    }
    //printf( "%d -> %d :n_send: %d \n", id, dest_id, n_send );
    
    reAllociBuffer( n_send[dest_id] );
        
    for( i = 0; i < n_send[dest_id]; ++i )
    {
        ibuffer[i] = mapping[i];
    }
        
    MPI_Send( ibuffer, n_send[dest_id], MPI_INT, dest_id, 99, MPI_COMM_WORLD );
           
}

void WallComm::recv_data( int source_id )
{

    //printf( "Starting Receive ( %d ) \n", id );

    MPI_Status status; 
    
    int i,j;
    //receive the number of particles 
    MPI_Recv( &n_recv[source_id], 1, MPI_INT, source_id, 99, MPI_COMM_WORLD, &status );
        
    if( n_recv[source_id] == 0 ) return;
    
    
    reAllociBuffer( n_recv[source_id] );
    MPI_Recv( ibuffer, n_recv[source_id] , MPI_INT, source_id, 99, MPI_COMM_WORLD, &status );
    
    for( i = 0; i < n_recv[source_id]; ++i )
    {
        local_ids.push_back( ibuffer[i] );
	cpu_ids.push_back( source_id );
    }
	
}

void WallComm::self_comm()
{
    
    std::vector<int>& mapping = elementCPUMapping_[id];
    
    for( int i = 0; i < mapping.size(); ++i )
    {
	local_ids.push_back( mapping[i] );
	cpu_ids.push_back( id );
    }
    
}

int WallComm::coord2bin(const double *x, int ic) const
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
      //if( xx[i] < bounds__[2*i] )   xx[i] = bounds__[2*i+1]-bounds__[2*i] + xx[i]; 
      //if( xx[i] > bounds__[2*i+1] ) xx[i] = bounds__[2*i]-bounds__[2*i+1] + xx[i];  

      iCell_lo = (xx[i] - this->coords_[6*ic+2*i]);
      iCell_hi = (xx[i] - this->coords_[6*ic+2*i+1]);
      
      if( iCell_lo >= -asciiPrecission && iCell_hi <= asciiPrecission )
	  iCell++;
  }

  return iCell;
  
}

vector<double> WallComm::parse_vector( string line, int i )
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
	        printf("Parsing failed due to unknown number!");
	    }
	    
	}else{
	    ++j;
	}
    }
    
    return parsed;
    
}

int WallComm::parse_int( string line )
{
   
    int j = 0;
    string nbuffer;
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

double WallComm::parse_double( string line, int* i, bool& success )
{
    int j = i[0];
    int index = 0;
    string nbuffer;
    
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

bool WallComm::is_digit( char c )
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

void WallComm::get_face_vector( vector<double*>& face, vector<double> parsed )
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










