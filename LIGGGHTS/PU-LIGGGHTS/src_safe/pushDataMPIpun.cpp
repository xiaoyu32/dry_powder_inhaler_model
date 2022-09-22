#include "cfd_datacoupling.h"
//#include "multisphere_parallel.h"
#include "error.h"
#include "properties.h"
#include "mpi.h"
#include "mpi_comm.h"
#include <vector>
#include <cstdlib>
#include "pushDataMPIpun.h"
#include "atom.h"
#include "cfd_datacoupling_mpi_pun.h"
#include "particleCommunicator.h"

#include  <iostream>

using namespace std;
using namespace LAMMPS_NS;

PushDataMPIpun::PushDataMPIpun( CfdDatacoupling * cfd_coupling_) : MPI_Communicator::MPI_Communicator()
{
    cfd_coupling = static_cast<CfdDatacouplingMPIpun*>(cfd_coupling_);
    cpu_mapping = cfd_coupling->cpu_mapping;
    init_send = false;
        
}

PushDataMPIpun::~PushDataMPIpun()
{
    if( buffer ) free( buffer );
    if( ibuffer ) free( ibuffer );
  
    buffer = NULL;
    ibuffer = NULL;
    init_send = false;
}

void PushDataMPIpun::do_comm(const char *name,const char *type,void *&to, const char* dataType)
{
    
    //printf( "Communicating: %s, %s  %d \n", name, type, id );
    
    fdata = cfd_coupling->find_push_property( name, type, n_dem, len );
    
    //get pointer to Push data to
    tdata = to; //TODO: change this to local one
    
    if( strncmp( dataType, "int", 3 ) == 0 ){
       //std::cout<<"int data"<<std::endl;
       data_type = MPI_INT;
    }else if( strncmp( dataType, "double", 6 ) == 0 ){
       //std::cout<<"Double data"<<std::endl;
       data_type = MPI_DOUBLE;
    }else
    {
       std::cerr<<"Unknown data type: "<<dataType<<std::endl;
       data_type = MPI_DOUBLE;
    }
       
    //std::cout<<" Name = "<<name
    //         <<" Type = "<<dataType
    //	     <<" len = "<<len<<std::endl;


    PushDataMPIpun::init_data_transfer();  
    
    //Push the data
    MPI_Communicator::do_bcast_comm();
    
}

void PushDataMPIpun::send_data( int dest_id )
{
    
    //init_send = true;    
    std::vector<int>& cpu_indexes = cfd_coupling->particleComm()->getCPUindexesDEM( dest_id );
        
    int n_send = cfd_coupling->particleComm()->nsend( dest_id );
    
    if( cpu_indexes.size() != n_send ) 
    {
        printf( "Error: inconsistent number detected in PushDataMPIpun::send_data( int dest_id )\n" );
    }
       
    int i,j;
     
    if( n_send == 0 )
    {
	return;
    }
    
    int n_bytes = n_send * len * sizeof( double ); 
    
    reAllocBuffer( n_bytes );
    
    if( data_type == MPI_DOUBLE ) // double data 
    {               
       double* buff2 = (double*)buffer;

       if( len == 1 ){

	  double* fbuff = (double*)fdata;

	  for( i = 0 ; i < n_send; ++i )
	  {
	      buff2[i] = fbuff[ cpu_indexes[i] ];   
	  }

       }else{

	  double** fbuff = (double**)fdata;

	  for( i = 0 ; i < n_send; ++i )
	  {
              for( j = 0; j < len; ++j )
	      {
		  buff2[i*len+j] = fbuff[ cpu_indexes[i] ][j];
	      }
	  }

       } 

       MPI_Send( buff2, n_send * len, data_type, dest_id, 99, MPI_COMM_WORLD );
       
    }else{ // int data
              
       int* buff2 = (int*)buffer;

       if( len == 1 ){

	  int* fbuff = (int*)fdata;
	
	  for( i = 0 ; i < n_send; ++i )
	  {
	      buff2[i] = fbuff[ cpu_indexes[i] ];   	      
	  }

       }else{

	  int** fbuff = (int**)fdata;

	  for( i = 0 ; i < n_send; ++i )
	  {
              for( j = 0; j < len; ++j )
	      {
		  buff2[i*len+j] = fbuff[ cpu_indexes[i] ][j];
	      }
	  }

       } 

       MPI_Send( buff2, n_send * len, data_type, dest_id, 99, MPI_COMM_WORLD );    
    }
       
}

void PushDataMPIpun::recv_data( int source_id )
{

    //printf( "Starting Receive ( %d ) \n", id );

    MPI_Status status; 
    
    int i,j;
        
    int n_recv = cfd_coupling->particleComm()->nrecv( source_id );	

    if( n_recv == 0 ) return;
    
    int n_bytes = n_recv* len * sizeof( double );
    
    reAllocBuffer( n_bytes );
    MPI_Recv( buffer, n_recv * len, data_type, source_id, 99, MPI_COMM_WORLD, &status );
    
    
    if( data_type == MPI_DOUBLE )
    {
       double ** tbuff = (double**)tdata;
       double * buff2 = (double*)buffer;
       
       if( len == 1 ){

	  for( i = 0 ; i < n_recv; ++i )
	  {
              //printf( "Index: %d \n", index );
	      tbuff[index][0] = buff2[i];
	      ++index;
	  }

       }else{

	  for( i = 0 ; i < n_recv; ++i )
	  {
              for( j = 0; j < len; ++j ){
		  tbuff[index][j] = buff2[i*len + j];
	      }
	      ++index;
	  }
       } 
    }else{ // int data
    
       int ** tbuff = (int**)tdata;
       int * buff2 = (int*)buffer;

       if( len == 1 ){

	  for( i = 0 ; i < n_recv; ++i )
	  {
              //printf( "Index: %d \n", index );
	      tbuff[index][0] = buff2[i];
	      ++index;
	  }

       }else{

	  for( i = 0 ; i < n_recv; ++i )
	  {
              for( j = 0; j < len; ++j ){
		  tbuff[index][j] = buff2[i*len + j];
	      }
	      ++index;
	  }
       }     
    }
     	
}

void PushDataMPIpun::self_comm()
{
    
    std::vector<int>& cpu_indexes = cfd_coupling->particleComm()->getCPUindexesDEM( id );
    
    if( data_type == MPI_DOUBLE )
    {
       double** tbuff = (double**)tdata;

       if( len == 1 ){

	  double* fbuff = (double*)fdata;

	  for( int i = 0; i < cpu_indexes.size(); ++i )
	  {
	      tbuff[index][0] = fbuff[ cpu_indexes.at(i) ];
	      ++index;
	  }
       }else{

	  double** fbuff = (double**)fdata;

	  for( int i = 0; i < cpu_indexes.size(); ++i )
	  {
	      for( int j = 0; j < len; ++j )
	      {
		 tbuff[index][j] = fbuff[ cpu_indexes.at(i) ][j];
	      }	
	      ++index;
	  }

       }
    }else{ // int data
       int** tbuff = (int**)tdata;

       if( len == 1 ){

	  int* fbuff = (int*)fdata;

	  for( int i = 0; i < cpu_indexes.size(); ++i )
	  {
	      tbuff[index][0] = fbuff[ cpu_indexes.at(i) ];	      
	      ++index;
	  }
       }else{

	  int** fbuff = (int**)fdata;

	  for( int i = 0; i < cpu_indexes.size(); ++i )
	  {
	      for( int j = 0; j < len; ++j )
	      {
		 tbuff[index][j] = fbuff[ cpu_indexes.at(i) ][j];
	      }	
	      ++index;
	  }

       }
    
    }
    
}

//find data that belongs to each CPU
void PushDataMPIpun::init_data_transfer()
{
    n_dem = cpu_mapping->n_dem;
    index = 0;
}


void PushDataMPIpun::end_of_comm()
{
     
     /*if( index != cpu_mapping->n_cfd )
     {
         printf( "Error! (6) \n" ); 
     }*/
     
}









