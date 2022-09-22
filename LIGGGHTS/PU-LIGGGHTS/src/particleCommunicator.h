#ifndef PARTICLE_COMMUNICATOR_PUN_H
#define PARTICLE_COMMUNICATOR_PUN_H


#include "mpi_comm.h"
#include "fix_eulerian_CPUs.h"
#include "cfd_datacoupling_mpi_pun.h"

namespace LAMMPS_NS {

/*
 * Establishes particle MPI task topology, i.e which DEM particle goes to which CFD task and vice versa
 * Used in the PUN cfd data coupling class
 * Jari Kolehmainen
 */

class particleCommunicator : public MPI_Communicator {

  public :

   particleCommunicator(class CfdDatacouplingMPIpun *);
   ~particleCommunicator();

   // Main communication
   void do_comm();
   
   inline int nsend(int i) const
   {
       return n_send[i];
   }
   
   inline int nrecv( int i ) const
   {
       return n_recv[i];
   }
   
   // -- returns list of particle (DEM) indices that go to CFD task i --
   inline std::vector<int>& getCPUindexesDEM( int i )
   {
       return cpu_indexes_dem[i];
   }
   
   // -- return list of particle (CFD) indices that go to DEM task i --
   inline std::vector<int>& getCPUindexesCFD( int i )
   {
       return cpu_indexes_cfd[i];
   }
   
   inline bool isUpdated() const
   {
       return updated;
   }  
   
   inline void setUpdated( bool updated_)
   {
       updated = updated_;
   }
   
   void printCommunicationPattern();
   
   std::vector<int> local_ids;
   std::vector<int> cpu_ids;
   
   protected:
     FixEulerianCPUs* cpu_mapping;	
	
     //virtual functions for handling data sending and receiving
     virtual void send_data( int );	// int: destination
     virtual void recv_data( int );    // int: source
     virtual void self_comm();         // communication in own rank
     virtual void end_of_comm();
     
     //number of particles send for a given processor
     int* n_send;
     int* n_recv;
     
     bool updated;
     
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
     
     //list particle global ids (needs to be structured in the same way as the from pointer
     std::vector<int>* cpu_indexes_dem;      // list of dem indices going to cfd processor id
     std::vector<int>* cpu_indexes_cfd;      // list of cfd indices going to dem processor id
         
};

}
#endif
