#include "mpi_comm.h"
#include "fix_eulerian_CPUs.h"
#include "cfd_datacoupling_mpi_pun.h"
#include "particleCommunicator.h"
#include <cstdlib>

using namespace LAMMPS_NS;

particleCommunicator::particleCommunicator( CfdDatacouplingMPIpun * cfd_coupling_ ) : 
MPI_Communicator::MPI_Communicator(),
cfd_coupling( NULL ),
updated( false )
{
    cpu_indexes_dem = new vector<int>[ntask];
    cpu_indexes_cfd = new vector<int>[ntask];
    
    n_send = (int*)malloc( sizeof(int) * ntask );
    n_recv = (int*)malloc( sizeof(int) * ntask );
    cfd_coupling = static_cast<CfdDatacouplingMPIpun*>(cfd_coupling_); 
    cpu_mapping = cfd_coupling->cpu_mapping;
}



particleCommunicator::~particleCommunicator()
{
    if( cpu_indexes_dem ) delete[] cpu_indexes_dem;
    if( cpu_indexes_cfd ) delete[] cpu_indexes_cfd;
    if( buffer ) free( buffer );
    if( ibuffer ) free( ibuffer );
    if( n_send ) free( n_send );
    if( n_recv ) free( n_recv );
    
    n_send = NULL;    
    n_recv = NULL;
    cpu_indexes_dem = NULL;
    cpu_indexes_cfd = NULL;
    buffer = NULL;
    ibuffer = NULL;
}



void particleCommunicator::do_comm()
{
    // -- communicator is upto data, no need to re-establish particle MPI task topology --
    if( updated ) return;
    
    // initialize communication
    particleCommunicator::init_data_transfer();  
    
    // update topology
    MPI_Communicator::do_bcast_comm();
    
    //this->printCommunicationPattern();
    
    updated = true;
}

void particleCommunicator::send_data( int dest_id )
{
    
    n_send[dest_id] = cpu_indexes_dem[dest_id].size();//cpu_mapping->localtagCFD()[dest_id];
    
    int n_bytes = n_send[dest_id] * len * sizeof( double );    
    int i,j;
    
    double * buff2;
    
    //send the number of particles that dest_id will be receiving
    MPI_Send( &n_send[dest_id], 1, MPI_INT, dest_id, 99, MPI_COMM_WORLD ); 
     
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
        ibuffer[i] = cpu_indexes_dem[dest_id][i];
    }
        
    MPI_Send( ibuffer, n_send[dest_id], MPI_INT, dest_id, 99, MPI_COMM_WORLD);
           
}

void particleCommunicator::recv_data( int source_id )
{

    //printf( "Starting Receive ( %d ) \n", id );

    MPI_Status status; 
    
    int i,j;
    //receive the number of particles 
    MPI_Recv( &n_recv[source_id], 1, MPI_INT, source_id, 99, MPI_COMM_WORLD, &status );
        
    if( n_recv[source_id] == 0 ) return;
    
    int n_bytes = n_recv[source_id] * len * sizeof( double );
    
    reAllociBuffer( n_recv[source_id] );
    MPI_Recv( ibuffer, n_recv[source_id] , MPI_INT, source_id, 99, MPI_COMM_WORLD, &status );
    
    for( i = 0; i < n_recv[source_id]; ++i )
    {
        local_ids.push_back( ibuffer[i] );
	cpu_ids.push_back( source_id );
    }
	
}

void particleCommunicator::self_comm()
{
    
    n_send[id] = cpu_indexes_dem[id].size();
    n_recv[id] = cpu_indexes_dem[id].size();
    
    for( int i = 0; i < cpu_indexes_dem[id].size(); ++i )
    {
	local_ids.push_back( cpu_indexes_dem[id][i] );
	cpu_ids.push_back( id );
    }
   
}

//find data that belongs to each CPU
void particleCommunicator::init_data_transfer()
{
    local_ids.clear();
    cpu_ids.clear();
    n_dem = cpu_mapping->n_dem;
    index = 0;

    for( int i = 0; i < ntask; ++i )
    {
	cpu_indexes_dem[i].clear();
    }

    for( int j = 0; j < n_dem; ++j )
    {
        int ic = cpu_mapping->cpu_indexes[j];
	if( ic >= 0 ){
	    cpu_indexes_dem[ ic ].push_back(j);
	}else{
	    //printf( "Warning: There are particles not been sent to any processor! \n" );
	}
    }
 
    //PushDataMPIpu::debug();
    
}

// -- reconstruc the inverse mapping -> which cfd domain particles go to which dem process --
void particleCommunicator::end_of_comm()
{
    
    for( int i = 0; i < ntask; ++i )
    {
	cpu_indexes_cfd[i].clear();
    }
      
    for( int j = 0; j < local_ids.size(); ++j )
    {
	cpu_indexes_cfd[cpu_ids[j]].push_back(j);
    }         
      
}

void particleCommunicator::printCommunicationPattern()
{
    
    for( int itask = 0; itask < ntask; ++itask )
    {
    
        if( id == itask )
	{
	    printf( "ID: %d  with %d particles \n", id, n_dem );
	    printf( "MPI Task    " );
	    
	    for( int jtask = 0; jtask < ntask; ++jtask )
	        printf( "  %d  ", jtask );
	    printf( "\n" );
	    
	    printf( "Send        " );
	    for( int jtask = 0; jtask < ntask; ++jtask )
	    {
	        printf( " %d ", n_send[jtask] );
	    }
	    printf( "\n" );
	    
	    printf( "Received    " );
	    for( int jtask = 0; jtask < ntask; ++jtask )
	    {
	        printf( " %d ", n_recv[jtask] );
	    }
	    printf( "\n" );	    
        }
	
	MPI_Barrier( MPI_COMM_WORLD );
	
    }
    
}

