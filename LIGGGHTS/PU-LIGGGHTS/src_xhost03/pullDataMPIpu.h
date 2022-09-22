#ifndef PULL_DATA_MPI_PU_H
#define PULL_DATA_MPI_PU_H

#include "cfd_datacoupling.h"
#include "error.h"
#include "properties.h"
#include "mpi.h"
#include "mpi_comm.h"
#include <vector>

#define INT_BUFFER_SIZE 20

namespace LAMMPS_NS {

class PullDataMPIpu : public MPI_Communicator {
 public:
  PullDataMPIpu(class CfdDatacoupling *);
  ~PullDataMPIpu();

   // Main communication
   virtual void do_comm(const char *name,const char *type,void *&from);
   
   std::vector<int> local_ids; //global ids corresponding to from pointers location
   std::vector<int> cpu_ids; //cpu ranks corresponding to from pointer location  
   
   protected:

     //virtual functions for handling data sending and receiving
     virtual void send_data( int );	// int: destination
     virtual void recv_data( int );    // int: source
     virtual void self_comm();         // communication in own rank
     
   private:
     void init_data_transfer();
     void reAllocIntBuffer(int);

     CfdDatacoupling * cfd_coupling;

     //local number of particles on the CFD side     
     int n_cfd;
      
     //local number of particles on the DEM side    
     int n_dem;

     //data vector length
     int len;
     
     //list particle global ids (needs to be structured in the same way as the from pointer  
     std::vector<int>* cpu_indexes;
     int * cpu_dem_ids; //global ids of the received particles	
     int cpu_dem_ids_length;
};
}
#endif
