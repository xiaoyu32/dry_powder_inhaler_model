/* ----------------------------------------------------------------------
   Jari Kolehmainen
------------------------------------------------------------------------- */

/*
   C or Fortran style library interface to LAMMPS
   new LAMMPS-specific functions can be added
*/


#ifndef LIGGGHTS_MPI_COMM_H
#define LIGGGHTS_MPI_COMM_H

#include "mpi.h"

#define BUFFER_CAPACITY 640000

#define IBUFFER_CAPACITY 100000

class MPI_Communicator{

   public:
 
      MPI_Communicator(); 
 
      ~MPI_Communicator();

      virtual void do_comm(const char *name,const char *type,void *&from );
    
   protected:
   
      //call send_data and recv_data between every mpitask
      void do_bcast_comm();
   
      //virtual functions for handling data sending and receiving
      virtual void send_data( int );	// int: destination
      virtual void recv_data( int );    // int: source
      virtual void self_comm();         // communication in own rank
      virtual void end_of_comm();	//end of communication
      
      //reallocate the communication buffer if necessary using bytes
      void reAllocBuffer(int);
	
      //reallocate int and double buffer for number of elements	
      void reAllociBuffer(int);	
	
      MPI_Datatype data_type;

      //task id
      int id;

      //number of mpi tasks
      int ntask;

      //data pointer
      void* fdata;
      void* tdata;
     
      void* buffer;
      int* ibuffer;
      
      int sbuffer; 
      int sibuffer;     

};


#endif 
