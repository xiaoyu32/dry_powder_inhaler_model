#ifndef PULL_DATA_MPI_PUN_H
#define PULL_DATA_MPI_PUN_H
#include "cfd_datacoupling.h"
#include "cfd_datacoupling_mpi_pun.h"
//#include "multisphere_parallel.h"
#include "error.h"
#include "properties.h"
#include "mpi.h"
#include "mpi_comm.h"
#include <vector>

#define INT_BUFFER_SIZE 20

namespace LAMMPS_NS {

class PullDataMPIpun : public MPI_Communicator 
{

 public:
  PullDataMPIpun(class CfdDatacoupling *);
  ~PullDataMPIpun();

   // Main communication
   virtual void do_comm(const char *name,const char *type,void *&from, const char* dataType);
   
   std::vector<int> local_ids; //global ids corresponding to from pointers location
   std::vector<int> cpu_ids; //cpu ranks corresponding to from pointer location  
   
   protected:

     //virtual functions for handling data sending and receiving
     virtual void send_data( int );	// int: destination
     virtual void recv_data( int );    // int: source
     virtual void self_comm();         // communication in own rank

     //status object for nonblocking communication
     MPI_Request destReqs;
     MPI_Request send_reqs[2];
     MPI_Status send_stats[2];
     
     
     // check if there has been isend calll
     bool init_send;
     
   private:
     void init_data_transfer();
     void reAllocIntBuffer(int);

     CfdDatacouplingMPIpun * cfd_coupling;

     //local number of particles on the CFD side     
     int n_cfd;
      
     //local number of particles on the DEM side    
     int n_dem;

     //data vector length
     int len;

};
}
#endif
