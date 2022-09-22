#ifndef PUSH_WALL_MPI_PU_H
#define PUSH_WALL_MPI_PU_H

#include "cfd_datacoupling.h"
#include "error.h"
#include "properties.h"
#include "mpi.h"
#include "mpi_comm.h"
#include "wallComm.h"
#include <vector>

namespace LAMMPS_NS {

class PushWallMPI : public MPI_Communicator {
 public:
  PushWallMPI( const class WallComm* wallComm_ );
  ~PushWallMPI();

   // Main communication
   void do_comm(void*& from, void *&to, int len_, int nvec_ );

   protected:
   
     const WallComm* wallComm;	
	
     //virtual functions for handling data sending and receiving
     virtual void send_data( int );	// int: destination
     virtual void recv_data( int );    // int: source
     virtual void self_comm();         // communication in own rank
     
     //status object for nonblocking communication
     MPI_Request destReqs;
     MPI_Request send_reqs;
     MPI_Status send_stats;
     
   private:
     
     //data vector length
     int len;
     
     //number of vectors
     int nvec;
     
     //index of CFD side particles
     int index;
     
};
}

#endif
