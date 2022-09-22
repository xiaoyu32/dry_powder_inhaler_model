#include "cfd_datacoupling.h"
#include "error.h"
#include "properties.h"
#include "mpi.h"
#include "mpi_comm.h"
#include <vector>
#include <cstdlib>
#include "pullDataMPIpu.h"
#include "atom.h"

using namespace std;
using namespace LAMMPS_NS;

PullDataMPIpu::PullDataMPIpu( CfdDatacoupling * cfd_coupling_) : MPI_Communicator::MPI_Communicator()
{
    cfd_coupling = cfd_coupling_;
    cpu_indexes = new vector<int>[ntask];
    cpu_dem_ids = (int*)malloc( INT_BUFFER_SIZE * sizeof(int) );
    cpu_dem_ids_length = INT_BUFFER_SIZE;
}

PullDataMPIpu::~PullDataMPIpu()
{
    if( cpu_indexes ) delete[] cpu_indexes;
    if( buffer ) free( buffer );
    if( cpu_dem_ids ) free( cpu_dem_ids );
    
    cpu_indexes = NULL;
    buffer = NULL;
    cpu_dem_ids = NULL;
    
}

void PullDataMPIpu::do_comm(const char *name,const char *type,void *&from)
{
        
    fdata = from;
    
    /*if( len == 3 )
    {
       double** ffdata = (double**)fdata;
       
       //for( int i =0; i < n_cfd; ++i )
       printf( "n_cfd = %d \n", n_cfd );
       printf( "From communication: %e %e %e \n", ffdata[0][0], ffdata[0][1], ffdata[0][2] );
    }*/
        
    //get pointer to pull data to
    tdata = cfd_coupling->find_pull_property( name, type, n_dem, len ); //TODO: change this to local one
        
    //TODO: do something smarter
    data_type = MPI_DOUBLE;

    init_data_transfer();    
    
    MPI_Barrier( MPI_COMM_WORLD );
    
    //pull the data
    MPI_Communicator::do_bcast_comm();
    
    MPI_Barrier( MPI_COMM_WORLD );
    
    /*if( len > 1 )
    {
        
        double** ttdata = (double**)tdata;
	printf( "After communication: %e %e %e \n", ttdata[0][0], ttdata[0][1], ttdata[0][2] );
    }*/
    
}

void PullDataMPIpu::send_data( int dest_id )
{

    int i_cfd = 0;
    int n_send = cpu_indexes[dest_id].size();
    int n_bytes = n_send * len * sizeof( double );    
    double** fbuff = (double**)fdata;
    
    int * buff;
    double * buff2;
    
    //send the number of particles that dest_id will be receiving
    MPI_Send( &n_send, 1, MPI_INT, dest_id, 99, MPI_COMM_WORLD );
    
    if( n_send == 0 ) return;
    
    //make sure the buffer size is adequote
    reAllocBuffer( n_bytes );    

    buff = (int*)buffer;
    
    //write particle global ids into the buffer
    for( int i = 0; i < n_send; ++i )
    {
	i_cfd = cpu_indexes[dest_id][i];
	buff[i] = local_ids[i_cfd];
    }
    
    //send buffer
    MPI_Send( buff, n_send, MPI_INT, dest_id, 99, MPI_COMM_WORLD );    
    
    buff2 = (double*)buffer;
    
    //write particle global ids into the buffer
    for( int i = 0; i < n_send; ++i )
    {
  	
	i_cfd = cpu_indexes[dest_id][i];
	
	for( int k = 0; k < len; ++k )
	{
	    buff2[k + i*len] = fbuff[i_cfd][k];
	}
    }
		
    //send buffer
    MPI_Send( buffer, n_send*len, data_type, dest_id, 99, MPI_COMM_WORLD );
    
}

void PullDataMPIpu::recv_data( int source_id )
{
    MPI_Status status; 
    int n_recv, i_local, m;

    //send the number of particles that dest_id will be receiving
    MPI_Recv( &n_recv, 1, MPI_INT, source_id, 99, MPI_COMM_WORLD, &status );

    if( n_recv == 0 ) return;

    int n_bytes = n_recv * len * sizeof( double );
    
    //make sure the buffer size is adequote
    reAllocIntBuffer( n_recv );
    reAllocBuffer( n_bytes );  
    
    MPI_Recv( cpu_dem_ids, n_recv, MPI_INT, source_id, 99, MPI_COMM_WORLD, &status );
    MPI_Recv( buffer, n_recv*len, data_type, source_id, 99, MPI_COMM_WORLD, &status );
    

    //insert particles to data pointer 
    double** tbuff = (double**)tdata;
    double* buff = (double*) buffer;    

    //read data from the buffer
    if( len == 1 )
    {

       double* tbuff = (double*)tdata;
       for( int i = 0; i < n_recv; ++i )
       {
	   i_local = cpu_dem_ids[i];
	   tbuff[i_local] = buff[i]; 
       }
    }else{
       double** tbuff = (double**)tdata;
       for( int i = 0; i < n_recv; ++i )
       {
	   i_local = cpu_dem_ids[i];
	   for( int k = 0; k < len; ++k )
	   {
	       tbuff[i_local][k] = buff[k + i*len]; 
    	   }
       }
    
    }

}


void PullDataMPIpu::self_comm()
{
    
    int j;
    int i_cfd;     
    double** fbuff = (double**)fdata;
        
    if( len == 1 )
    {

       double* tbuff = (double*)tdata;
       for( int i = 0; i < cpu_indexes[id].size(); ++i )
       {
           i_cfd = cpu_indexes[id][i];
	   j = local_ids[i_cfd];

	   tbuff[j] = fbuff[ i_cfd ][0]; //fbuff[ i_cfd ]; //

       }
    
    }else
    { 
              
       double** tbuff = (double**)tdata;
       
       for( int i = 0; i < cpu_indexes[id].size(); ++i )
       {
           i_cfd = cpu_indexes[id][i];
	   j = local_ids[i_cfd];
	   	   
	   for( int k = 0; k < len; ++k )
	   {
	       tbuff[j][k] = fbuff[ i_cfd ][k]; //fbuff[ k+i_cfd*len ];
	   }
	   
	   /*if( len > 1 && j == 0 )
	   {
	    
	       printf( "Inside communication: %e %e %e \n", tbuff[j][0], tbuff[j][1], tbuff[j][2] );
 	   }*/
       }
    
    }

}

//find data that belongs to each CPU
void PullDataMPIpu::init_data_transfer()
{
    
    for( int i = 0; i < ntask; ++i )
    {
	cpu_indexes[i].clear();
    }

    for( int j = 0; j < local_ids.size(); ++j )
    {
	cpu_indexes[cpu_ids[j]].push_back(j);
    }

}

void PullDataMPIpu::reAllocIntBuffer( int length )
{
    
    if( length > (cpu_dem_ids_length-1) ){
	cpu_dem_ids = (int*)realloc( cpu_dem_ids, sizeof(int)*(length) );
	cpu_dem_ids_length = length;
    }

}






