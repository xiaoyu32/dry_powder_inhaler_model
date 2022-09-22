#ifndef PUSH_DATA_MPI_PUN_H
#define PUSH_DATA_MPI_PUN_H

#include "cfd_datacoupling.h"
//#include "multisphere_parallel.h"
#include "error.h"
#include "properties.h"
#include "mpi.h"
#include "mpi_comm.h"
#include <vector>
#include "fix_eulerian_CPUs.h"
#include "cfd_datacoupling_mpi_pu.h"

#define INT_BUFFER_SIZE 20

namespace LAMMPS_NS {

class PushDataMPIpun : public MPI_Communicator {
 public:
  PushDataMPIpun(class CfdDatacoupling *);
  ~PushDataMPIpun();

   // Main communication
   virtual void do_comm(const char *name,const char *type,void *&to, const char* dataType );
   
   protected:
     FixEulerianCPUs* cpu_mapping;	
	
     //virtual functions for handling data sending and receiving
     virtual void send_data( int );	// int: destination
     virtual void recv_data( int );    // int: source
     virtual void self_comm();         // communication in own rank
     virtual void end_of_comm();
          

     //status object for nonblocking communication
     MPI_Request destReqs;
     MPI_Request send_reqs[2];
     MPI_Status send_stats[2];
     
     // check if there has been isend calll
     bool init_send;
     
     
   private:
   
     void init_data_transfer();

     class CfdDatacouplingMPIpun * cfd_coupling;

     //local number of particles on the CFD side     
     int n_cfd;
      
     //local number of particles on the DEM side    
     int n_dem;

     //data vector length
     int len;
     
     //index of CFD side particles
     int index;
     

};
}
#endif
