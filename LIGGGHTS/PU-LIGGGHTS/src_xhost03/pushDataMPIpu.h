#ifndef PUSH_DATA_MPI_PU_H
#define PUSH_DATA_MPI_PU_H

#include "cfd_datacoupling.h"
#include "error.h"
#include "properties.h"
#include "mpi.h"
#include "mpi_comm.h"
#include <vector>
#include "fix_eulerian_CPUs.h"
#include "cfd_datacoupling_mpi_pu.h"

#define INT_BUFFER_SIZE 20

namespace LAMMPS_NS {

class PushDataMPIpu : public MPI_Communicator {
 public:
  PushDataMPIpu(class CfdDatacoupling *);
  ~PushDataMPIpu();

   // Main communication
   virtual void do_comm(const char *name,const char *type,void *&to);

   protected:
     FixEulerianCPUs* cpu_mapping;	
	
     //virtual functions for handling data sending and receiving
     virtual void send_data( int );	// int: destination
     virtual void recv_data( int );    // int: source
     virtual void self_comm();         // communication in own rank
     virtual void end_of_comm();
     
     void debug();
     
   private:
     void init_data_transfer();

     class CfdDatacouplingMPIpu * cfd_coupling;

     //local number of particles on the CFD side     
     int n_cfd;
      
     //local number of particles on the DEM side    
     int n_dem;

     //data vector length
     int len;
     
     //index of CFD side particles
     int index;
     
     //list particle global ids (needs to be structured in the same way as the from pointer
     std::vector<int>* cpu_indexes;

};
}
#endif
