#include "pushWallMPI.h"
#include "wallComm.h"

using namespace LAMMPS_NS;

PushWallMPI::PushWallMPI( const WallComm* wallComm_ ) :
MPI_Communicator(),
wallComm( wallComm_ ),
len( 0 ),
nvec( 0 ),
index( 0 )
{}

PushWallMPI::~PushWallMPI()
{}

void PushWallMPI::do_comm( void*& from, void *&to, int len_, int nvec_ )
{
    
    len = len_;
    nvec = nvec_;
    
    fdata = from;
    tdata = to;
    
    data_type = MPI_DOUBLE;
    
    index = 0;
    
    MPI_Communicator::do_bcast_comm();
    
}

void PushWallMPI::send_data( int dest_id )
{
    
    //init_send = true;    
    const std::vector<int>& elementMap = wallComm->elementCPUMapping( dest_id );
    
    //cpu_mapping->localtagCFD()[dest_id];
    
    int n_bytes = wallComm->nsend(dest_id) * len * nvec * sizeof( double );    
    int i,j,k;
    
    double * buff2;
    
    if( wallComm->nsend(dest_id) == 0 )
    {
	return;
    }

    reAllocBuffer( n_bytes );
                   
    buff2 = (double*)buffer;

    double*** fbuff = (double***)fdata;

    for( i = 0 ; i < wallComm->nsend(dest_id); ++i )
    {
	for( j = 0; j < len; ++j )
	{
	    for( int k = 0; k < nvec; ++k )
	    {
	       buff2[i*len*nvec + j*nvec + k] = fbuff[ elementMap[i] ][j][k];
	    }
	}
    }
    
    
    
    MPI_Isend( buff2, wallComm->nsend(dest_id) * len * nvec, MPI_DOUBLE, dest_id, 99, MPI_COMM_WORLD, &send_reqs );

    MPI_Waitall( 1, &send_reqs, &send_stats );
       
}

void PushWallMPI::recv_data( int source_id )
{

    //printf( "Starting Receive ( %d ) \n", id );

    MPI_Status status; 
    
    int i,j,k;
        
    if( wallComm->nrecv(source_id) == 0 ) return;
    
    double *** tbuff = (double***)tdata;
    int n_bytes = wallComm->nrecv(source_id) * len * nvec * sizeof( double );
    
    reAllocBuffer( n_bytes );
    MPI_Recv( buffer, wallComm->nrecv(source_id) * len * nvec, MPI_DOUBLE, source_id, 99, MPI_COMM_WORLD, &status );
      
    double * buff2 = (double*)buffer;

    for( i = 0 ; i < wallComm->nrecv(source_id); ++i )
    {
        for( j = 0; j < len; ++j )
	{
	    for( k = 0; k < nvec; ++k )
	    {
	       tbuff[index][j][k] = buff2[i*len*nvec + j*nvec + k];
	    }
	}
	++index;
    }
    
}

void PushWallMPI::self_comm()
{
    
    const std::vector<int>& elementMap = wallComm->elementCPUMapping( id );
    
    double*** tbuff = (double***)tdata;
    double*** fbuff = (double***)fdata;
    
    for( int i = 0; i < elementMap.size(); ++i )
    {
	for( int j = 0; j < len; ++j )
	{
	   for( int k = 0; k < nvec; ++k )
	       tbuff[index][j][k] = fbuff[ elementMap[i] ][j][k];
	}	
	++index;
    }
    
}



