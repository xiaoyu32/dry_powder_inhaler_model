#include "cfd_datacoupling.h"
//#include "multisphere_parallel.h"
#include "error.h"
#include "properties.h"
#include "mpi.h"
#include "mpi_comm.h"
#include <vector>
#include <cstdlib>
#include "pushDataMPIpu.h"
#include "atom.h"
#include "cfd_datacoupling_mpi_pu.h"

using namespace std;
using namespace LAMMPS_NS;

PushDataMPIpu::PushDataMPIpu( CfdDatacoupling * cfd_coupling_) : MPI_Communicator::MPI_Communicator()
{
    cfd_coupling = static_cast<CfdDatacouplingMPIpu*>(cfd_coupling_);
    cpu_mapping = cfd_coupling->cpu_mapping;
    cpu_indexes = new vector<int>[ntask];
}

PushDataMPIpu::~PushDataMPIpu()
{
    if( cpu_indexes ) delete[] cpu_indexes;
    if( buffer ) free( buffer );
    
    cpu_indexes = NULL;
    buffer = NULL;
}

void PushDataMPIpu::do_comm(const char *name,const char *type,void *&to)
{
    
    //printf( "Communicating: %s, %s  %d \n", name, type, id );
    
    fdata = cfd_coupling->find_push_property( name, type, n_dem, len );
    
    //get pointer to Push data to
    tdata = to; //TODO: change this to local one
    
    //TODO: do something smarter
    data_type = MPI_DOUBLE;


    PushDataMPIpu::init_data_transfer();  
    
    //Push the data
    MPI_Communicator::do_bcast_comm();
    
}

void PushDataMPIpu::send_data( int dest_id )
{
    
    int n_send = cpu_indexes[dest_id].size();//cpu_mapping->localtagCFD()[dest_id];
    
    int n_bytes = n_send * len * sizeof( double );    
    int i,j;
    
    int* buff;
    double* buff2;
    
    //send the number of particles that dest_id will be receiving
    MPI_Send( &n_send, 1, MPI_INT, dest_id, 99, MPI_COMM_WORLD );
    
     
    if( n_send == 0 ) return;
    
    //printf( "%d -> %d :n_send: %d \n", id, dest_id, n_send );
    
    reAllocBuffer( n_bytes );
    
    buff = (int*) buffer;
    for( i = 0; i < n_send; ++i )
    {
        buff[i] = cpu_indexes[dest_id][i];
    }
        
    MPI_Send( buff, n_send, MPI_INT, dest_id, 99, MPI_COMM_WORLD );
                   
    buff2 = (double*)buffer;
    

    if( len == 1 ){
    
       double* fbuff = (double*)fdata;

       for( i = 0 ; i < n_send; ++i )
       {
	   buff2[i] = fbuff[ cpu_indexes[dest_id][i] ];   

       }
    }else{

       double** fbuff = (double**)fdata;
    
       for( i = 0 ; i < n_send; ++i )
       {
           for( j = 0; j < len; ++j ){
	       buff2[i*len+j] = fbuff[ cpu_indexes[dest_id][i] ][j];
	   }
       }
    
    } 
    
    MPI_Send( buff2, n_send*len, MPI_DOUBLE, dest_id, 99, MPI_COMM_WORLD );
        
}

void PushDataMPIpu::recv_data( int source_id )
{

    //printf( "Starting Receive ( %d ) \n", id );

    MPI_Status status; 
    int n_recv, i_global, m;
    int i,j;
    //receive the number of particles 
    MPI_Recv( &n_recv, 1, MPI_INT, source_id, 99, MPI_COMM_WORLD, &status );
        
    if( n_recv == 0 ) return;
    
    double ** tbuff = (double**)tdata;
    int n_bytes = n_recv * len * sizeof( double );
    reAllocBuffer( n_bytes );
    
    MPI_Recv( buffer, n_recv , MPI_INT, source_id, 99, MPI_COMM_WORLD, &status );
    
    int* buff = (int*)buffer;
    for( i = 0; i < n_recv; ++i )
    {
        cfd_coupling->pull_->local_ids.push_back( buff[i] );
	cfd_coupling->pull_->cpu_ids.push_back( source_id );
    }
    
    
    MPI_Recv( buffer, n_recv * len, MPI_DOUBLE, source_id, 99, MPI_COMM_WORLD, &status );
    
   double * buff2 = (double*)buffer;

   //printf( "%d -> %d :n_recv: %d \n", source_id, id, n_recv );
   
   //printf( "Writing Data... ( %d ) \n", id );
   
   //sanity check
   //if( index + n_recv > cpu_mapping->n_cfd ){
   //    printf( "Error! (5) \n" );
       //cpu_mapping->error->all(FLERR,"TOO MANY PARTICLES!!!");
   //}
   
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
        
    //printf( "Data Writen ( %d ) \n", id );
	
}

void PushDataMPIpu::self_comm()
{
    
    
    double** tbuff = (double**)tdata;

    if( len == 1 ){
       
       double* fbuff = (double*)fdata;

       for( int i = 0; i < cpu_indexes[id].size(); ++i )
       {
	   tbuff[index][0] = fbuff[ cpu_indexes[id].at(i) ];
	   //printf( "local index = %d \n",cpu_indexes[id].at(i) );
	   cfd_coupling->pull_->local_ids.push_back( cpu_indexes[id][i] );
	   cfd_coupling->pull_->cpu_ids.push_back( id );
	   ++index;
       }
    }else{

       double** fbuff = (double**)fdata;
    
       for( int i = 0; i < cpu_indexes[id].size(); ++i )
       {
	   for( int j = 0; j < len; ++j )
	   {
	      tbuff[index][j] = fbuff[ cpu_indexes[id].at(i) ][j];
	      //printf( "local index = %d \n",cpu_indexes[id].at(i) );
	      cfd_coupling->pull_->local_ids.push_back( cpu_indexes[id][i] );
	      cfd_coupling->pull_->cpu_ids.push_back( id );
	   }	
	   ++index;
       }
    
    }
    
}

//find data that belongs to each CPU
void PushDataMPIpu::init_data_transfer()
{
    cfd_coupling->pull_->local_ids.clear();
    cfd_coupling->pull_->cpu_ids.clear();
    n_dem = cpu_mapping->n_dem;
    index = 0;

    for( int i = 0; i < ntask; ++i )
    {
	cpu_indexes[i].clear();
    }

    for( int j = 0; j < n_dem; ++j )
    {
        int ic = cpu_mapping->cpu_indexes[j];
	if( ic >= 0 ){
	    cpu_indexes[ cpu_mapping->cpu_indexes[j] ].push_back(j);
	}else{
	    //printf( "Warning: There are particles not been sent to any processor! \n" );
	}
    }
 
    //PushDataMPIpu::debug();
    
}


void PushDataMPIpu::debug()
{
     
     int my_task[ntask];
     int n_dem_[ntask];
     
     for( int i = 0; i < ntask; ++i )
     {
         my_task[i] = 0;
	 n_dem_[i] = 0;
     }
     
     for( int i = 0;i < ntask; ++i )
     {
         my_task[i] = cpu_indexes[i].size();
     }	
     	
	
     MPI_Allreduce( my_task, n_dem_, ntask, MPI_INT, MPI_SUM, MPI_COMM_WORLD );	
	
     /*if( n_dem_[id] != cpu_mapping->n_cfd )
     {
         printf( "Error! (3) \n" );
         //cpu_mapping->error->all(FLERR,"Wrong number of Particles!!!!");
     }
     
     printf( "Current=%d \n", n_dem_[id] );
     printf( "n_cfd=%d \n", cpu_mapping->n_cfd );
	 
     for( int i = 0; i < ntask; ++i )
     {
	     
         printf( "ID( %d->%d ) : %d , %d \n ", id, i, 
 		cpu_indexes[i].size(), cpu_mapping->localtagCFD()[i] );
	     
     }*/
	 
     
	
}

void PushDataMPIpu::end_of_comm()
{
     
     /*if( index != cpu_mapping->n_cfd )
     {
         printf( "Error! (6) \n" ); 
     }*/
     
}









