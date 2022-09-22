/* ----------------------------------------------------------------------
 Jari Kolehmainen
------------------------------------------------------------------------- */

#include "mpi.h"
#include "mpi_comm.h"
#include <cstdlib>
#include <cstdio>

MPI_Communicator::MPI_Communicator()
{
    MPI_Comm_rank( MPI_COMM_WORLD, &id );
    MPI_Comm_size( MPI_COMM_WORLD, &ntask );
    
    fdata = NULL;
    tdata = NULL;

    buffer = (void*) malloc( BUFFER_CAPACITY );
    ibuffer = (int*) malloc( IBUFFER_CAPACITY * sizeof( int ) );
    sbuffer = BUFFER_CAPACITY;
    sibuffer = IBUFFER_CAPACITY;
}

MPI_Communicator::~MPI_Communicator()
{
   if( buffer ) free(buffer);
   if( ibuffer ) free( ibuffer );
   
   ibuffer = NULL;
   buffer = NULL;
}

void MPI_Communicator::reAllocBuffer( int size )
{
    if( size > sbuffer )
    {
        printf( "Warning: Reallocated MPI communication buffer! \n" );
	buffer = realloc( buffer, 2*size  );
	sbuffer = 2*size;
    }
}

void MPI_Communicator::reAllociBuffer( int size )
{
    if( size > sibuffer )
    {
        printf( "Warning: Reallocated MPI communication buffer! \n" );
	
 	ibuffer = (int*)realloc( ibuffer, 2*size * sizeof( int ) );
	sibuffer = 2*size;
	
    }
}


void MPI_Communicator::do_comm(const char *name,const char *type,void *&from){}

void MPI_Communicator::send_data( int dest_id ){}

void MPI_Communicator::self_comm(){}

void MPI_Communicator::recv_data( int source_id ){}

void MPI_Communicator::end_of_comm(){}

void MPI_Communicator::do_bcast_comm()
{
    
    int i, sid, did, ID, turn;
    
    //printf( "Bcasting Comm... %d \n", id );
    //MPI_Barrier( MPI_COMM_WORLD );	

    //handle the self communication
    self_comm();
    
     //loop over the mpi tasks
    for( i = 1; i <= ntask/2; ++i )
    {
       	//this is always 0 or 1
	ID = (id-id%i)/i;
	
	for( turn = 0; turn < 2; ++ turn )
	{
 	    
	    if( ID%2 == turn )
	    {
		
		sid = (id+i)%ntask;		
		send_data( sid );

		if( i < ntask/2 ) recv_data( sid );

	    }else {
		
		did = (ntask+id-i)%ntask;
		recv_data( did );
		
		if( i < ntask/2 ) send_data( did );

	    }
	    	
	}
	
    } 
    
    end_of_comm();
    
    //printf( "Bcasted Comm %d \n", id );
    //MPI_Barrier( MPI_COMM_WORLD );
    
}















