#include "pullWallMPI.h"
#include "wallComm.h"

using namespace LAMMPS_NS;

PullWallMPI::PullWallMPI( const WallComm* wallComm_ ) :
MPI_Communicator(),
wallComm( wallComm_ ),
len( 0 ),
nvec( 0 ),
index( 0 )
{}

PullWallMPI::~PullWallMPI()
{}

void PullWallMPI::do_comm( void*& from, void *&to, int len_, int nvec_ )
{
    
    len = len_;
    nvec = nvec_;
    
    fdata = from;
    tdata = to;
    
    data_type = MPI_DOUBLE;
    
    index = 0;
    
    MPI_Communicator::do_bcast_comm();
    
}

void PullWallMPI::send_data( int dest_id )
{
    
    //init_send = true;

    const std::vector<int>& elementMap = wallComm->elementDEMMapping( dest_id );

    int i_cfd = 0;

    if( wallComm->nrecv(dest_id) == 0 )
    {
	return;
    }

    int n_bytes = wallComm->nrecv(dest_id) * len * nvec * sizeof( double );    
    double*** fbuff = (double***)fdata;
    
    double * buff2;
            
    reAllocBuffer( n_bytes );
    
    buff2 = (double*)buffer;
    
    //write particle global ids into the buffer
    for( int i = 0; i <  wallComm->nrecv(dest_id); ++i )
    {
  	
	i_cfd = elementMap[i];
	
	for( int j = 0; j < len; ++j )
	    for( int k = 0; k < nvec; ++k )
	       buff2[k + j*nvec+ i*len*nvec] = fbuff[i_cfd][j][k];
	
    }
		
    // -- send buffer
    MPI_Send( buffer, wallComm->nrecv(dest_id)*len*nvec, data_type, dest_id, 99, MPI_COMM_WORLD );

}

void PullWallMPI::recv_data( int source_id )
{
    MPI_Status status; 
    int i_local, m;
    
    const std::vector<int>& elementMap = wallComm->elementCPUMapping( source_id );
    	
    if( wallComm->nsend(source_id) == 0 ) return;
    
    //printf( "Inside communication: %d \n", n_recv );
    
    int n_bytes =  wallComm->nsend(source_id) * len * nvec * sizeof( double );
    

    reAllocBuffer( n_bytes );
    
    double* buff = (double*) buffer;
    MPI_Recv( buff,  wallComm->nsend(source_id) * len * nvec, data_type, source_id, 99, MPI_COMM_WORLD, &status );
    

    //insert particles to data pointer 
    double*** tbuff = (double***)tdata;

    for( int i = 0; i < wallComm->nsend(source_id); ++i )
    {
	i_local = elementMap[i];
	
	for( int j = 0; j < len; ++j )
	   for( int k = 0; k < nvec; ++k )
	   {
	       tbuff[i_local][j][k] = buff[k + j * nvec + i*len*nvec]; 
    	   }
    }
    
}

void PullWallMPI::self_comm()
{
    
    int i_cfd,i_dem;     

              
    double*** tbuff = (double***)tdata;
    double*** fbuff = (double***)fdata;
    
    // indices of Dem map that go to processor 'id'
    const std::vector<int>& elementMapDem = wallComm->elementCPUMapping( id );
    
    // indices of Cfd map that go to processor 'id' 
    const std::vector<int>& elementMapCfd = wallComm->elementDEMMapping( id );
    
    if( elementMapDem.size() != elementMapCfd.size() )
    {
        printf( "Error: element maps are inconsistent! \n" );
    }
       
    for( int i = 0; i < elementMapDem.size(); ++i )
    {
	//i_cfd = elementMap[i];
	//i_dem = local_ids[i_cfd];
	i_dem = elementMapDem[i];
	i_cfd = elementMapCfd[i];
	
	for( int j = 0; j < len; ++j )
	   for( int k = 0; k < nvec; ++k )
	   {
	       tbuff[i_dem][j][k] = fbuff[ i_cfd ][j][k];
    	   }

    }

}













