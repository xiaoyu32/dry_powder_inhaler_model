#include "cfd_datacoupling.h"
#include "cfd_datacoupling_mpi_pun.h"
//#include "multisphere_parallel.h"
#include "error.h"
#include "properties.h"
#include "mpi.h"
#include "mpi_comm.h"
#include <vector>
#include <cstdlib>
#include "pullDataMPIpun.h"
#include "atom.h"

#include  <iostream>

using namespace std;
using namespace LAMMPS_NS;

PullDataMPIpun::PullDataMPIpun( CfdDatacoupling * cfd_coupling_) : MPI_Communicator::MPI_Communicator()
{
    cfd_coupling = static_cast<CfdDatacouplingMPIpun*>(cfd_coupling_);
    init_send = false;  
}

PullDataMPIpun::~PullDataMPIpun()
{
    if( buffer ) free( buffer );
    if( ibuffer ) free( ibuffer );
  
    buffer = NULL;
    ibuffer = NULL;
    init_send = false;
    
}

void PullDataMPIpun::do_comm(const char *name,const char *type,void *&from, const char* dataType)
{
        
    fdata = from;

    //get pointer to pull data to
    tdata = cfd_coupling->find_pull_property( name, type, n_dem, len ); //TODO: change this to local one
        
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

    init_data_transfer();    
        
    //pull the data
    MPI_Communicator::do_bcast_comm();
    
}

void PullDataMPIpun::send_data( int dest_id )
{
    
    //init_send = true;

    int i_cfd = 0;

    // nsend equals now to the items received from DEM side 
    int n_send = cfd_coupling->particleComm()->nrecv( dest_id );
    
    if( n_send == 0 )
    {
	return;
    }
    
    std::vector<int>& cpu_indexes = cfd_coupling->particleComm()->getCPUindexesCFD( dest_id );
 
     if( n_send != cpu_indexes.size() )
    {
        printf( "Error: inconsistent size detected in PullDataMPIpun::send_data( int dest_id )!\n" );
    }
    
    int n_bytes = n_send * len * sizeof( double );     
     
    if( data_type == MPI_DOUBLE )
    {   
       double** fbuff = (double**)fdata;

       double* buff2;

       reAllocBuffer( n_bytes );

       buff2 = (double*)buffer;

       //write particle global ids into the buffer
       for( int i = 0; i < n_send; ++i )
       {

	   i_cfd = cpu_indexes[i];

	   for( int k = 0; k < len; ++k )
	   {
	       buff2[k + i*len] = fbuff[i_cfd][k];
	   }
       }

       //send buffer
       MPI_Send( buffer, n_send*len, data_type, dest_id, 99, MPI_COMM_WORLD );
       
    }else{
       int** fbuff = (int**)fdata;

       int* buff2;

       reAllocBuffer( n_bytes );

       buff2 = (int*)buffer;

       //write particle global ids into the buffer
       for( int i = 0; i < n_send; ++i )
       {

	   i_cfd = cpu_indexes[i];

	   for( int k = 0; k < len; ++k )
	   {
	       buff2[k + i*len] = fbuff[i_cfd][k];
	   }
       }

       //send buffer
       MPI_Send( buffer, n_send*len, data_type, dest_id, 99, MPI_COMM_WORLD );
    
    }
}

void PullDataMPIpun::recv_data( int source_id )
{
    MPI_Status status; 
    int i_local, m;

    // -- receive the number of particles from cfd that was sent there --
    int n_recv = cfd_coupling->particleComm()->nsend( source_id );
     	
    if( n_recv == 0 ) return;
        
    int n_bytes = n_recv * len * sizeof( double );
    
    //make sure the buffer size is adequote
    
     
    reAllocBuffer( n_bytes );
    
    // -- put the particles in the same place where they were sent --
    std::vector<int>& cpu_indexes = cfd_coupling->particleComm()->getCPUindexesDEM( source_id );
    
    if( n_recv != cpu_indexes.size() )
    {
        printf( "Error: inconsistent size detected in PullDataMPIpun::recv_data( int source_id )!\n" );
    }    
    
    MPI_Recv( buffer, n_recv * len, data_type, source_id, 99, MPI_COMM_WORLD, &status );
    
    if( data_type == MPI_DOUBLE )
    {
	double* buff = (double*) buffer;

	//read data from the buffer
	if( len == 1 )
	{
	   double* tbuff = (double*)tdata;
	   for( int i = 0; i < n_recv; ++i )
	   {
	       i_local = cpu_indexes[i];
	       tbuff[i_local] = buff[i]; 
	   }
	}else{
	   double** tbuff = (double**)tdata;
	   for( int i = 0; i < n_recv; ++i )
	   {
	       i_local = cpu_indexes[i];
	       for( int k = 0; k < len; ++k )
	       {
		   tbuff[i_local][k] = buff[k + i*len]; 
    	       }
	   }

	}
    }else{
	int* buff = (int*) buffer;

	//read data from the buffer
	if( len == 1 )
	{
	   int* tbuff = (int*)tdata;
	   for( int i = 0; i < n_recv; ++i )
	   {
	       i_local = cpu_indexes[i];
	       tbuff[i_local] = buff[i]; 
	   }
	}else{
	   int** tbuff = (int**)tdata;
	   for( int i = 0; i < n_recv; ++i )
	   {
	       i_local = cpu_indexes[i];
	       for( int k = 0; k < len; ++k )
	       {
		   tbuff[i_local][k] = buff[k + i*len]; 
    	       }
	   }

	}
    
    }

}


void PullDataMPIpun::self_comm()
{
    
    int i_cfd, i_dem;     
    
    std::vector<int>& cpu_indexes_cfd = cfd_coupling->particleComm()->getCPUindexesCFD( id );
    std::vector<int>& cpu_indexes_dem = cfd_coupling->particleComm()->getCPUindexesDEM( id );
    
    if( cpu_indexes_cfd.size() != cpu_indexes_dem.size() ) 
    {
        printf( "Error: inconsistent sizes detected in PullDataMPIpun::self_comm()!\n" );
    }
    if( data_type == MPI_DOUBLE )
    {
       double** fbuff = (double**)fdata;

       if( len == 1 )
       {

	  double* tbuff = (double*)tdata;
	  for( int i = 0; i < cpu_indexes_cfd.size(); ++i )
	  {
              i_cfd = cpu_indexes_cfd[i];
	      i_dem = cpu_indexes_dem[i];

	      tbuff[i_dem] = fbuff[ i_cfd ][0]; //fbuff[ i_cfd ]; //

	  }

       }else
       { 

	  double** tbuff = (double**)tdata;

	  for( int i = 0; i < cpu_indexes_cfd.size(); ++i )
	  {
              i_cfd = cpu_indexes_cfd[i];
	      i_dem = cpu_indexes_dem[i];

	      for( int k = 0; k < len; ++k )
	      {
		  tbuff[i_dem][k] = fbuff[ i_cfd ][k]; //fbuff[ k+i_cfd*len ];
	      }

	  }

       }
    }else{
       int** fbuff = (int**)fdata;

       if( len == 1 )
       {

	  int* tbuff = (int*)tdata;
	  for( int i = 0; i < cpu_indexes_cfd.size(); ++i )
	  {
              i_cfd = cpu_indexes_cfd[i];
	      i_dem = cpu_indexes_dem[i];

	      tbuff[i_dem] = fbuff[ i_cfd ][0]; //fbuff[ i_cfd ]; //

	  }

       }else
       { 

	  int** tbuff = (int**)tdata;

	  for( int i = 0; i < cpu_indexes_cfd.size(); ++i )
	  {
              i_cfd = cpu_indexes_cfd[i];
	      i_dem = cpu_indexes_dem[i];

	      for( int k = 0; k < len; ++k )
	      {
		  tbuff[i_dem][k] = fbuff[ i_cfd ][k]; //fbuff[ k+i_cfd*len ];
	      }

	  }

       }
    
    }

}

//find data that belongs to each CPU
void PullDataMPIpun::init_data_transfer()
{}

void PullDataMPIpun::reAllocIntBuffer( int length )
{}






