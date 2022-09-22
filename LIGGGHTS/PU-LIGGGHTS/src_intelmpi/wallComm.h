/*
 * Jari Kolehmainen, 2019 
 * Class for handling MPI topology (which DEM wall elements should go to whcih processor) 
 * for wall communication
 */

#ifndef LMP_WALL_COMM_H
#define LMP_WALL_COMM_H

#include "mpi.h"
#include "mpi_comm.h"
#include <vector>
#include "fix_eulerian_CPUs.h"
#include "wallElementFinder.h"

// -- book keeping which wall element belongs to which DEM and CFD processor --

namespace LAMMPS_NS
{

class WallComm : public MPI_Communicator
{
    
    public :
    
    WallComm();
        
    ~WallComm();
    
    // -- dummy --
    virtual void do_comm(const char *name,const char *type,void *&from){}
    
    // initialize the CFD assignments (need to be called if DEM processors are changed
    
    void init( double** arr_, int nlocal );
    
    inline const int& ncfd() const
    {
        return this->n_cfd;
    }

    inline const int& ndem() const
    {
        return n_dem;
    }
    
    inline int ntasks() const
    {
       return ntask;
    }
    
    // return element map for processor id
    inline const std::vector<int>& elementCPUMapping( int id ) const
    {
       return elementCPUMapping_[id];
    }

    // return element map for processor id
    inline const std::vector<int>& elementDEMMapping( int id ) const
    {
       return elementDEMMapping_[id];
    }
    
    // number of elements sent from DEM to CFD
    inline int nsend( int id ) const
    {
	return n_send[id];
    }
    
    // number of elements received from DEM to CFD
    inline int nrecv( int id ) const
    {  
        return n_recv[id];
    } 
    
    inline const std::vector<int>& localIDS() const
    {
        return local_ids;
    }
    
    inline bool isActive() const
    {
        return active;
    }
    
    inline void setActive( bool active_ )
    {
        active = active_;
    }
    
    
    protected :
    
    WallElementFinder** wallElementFinder_;
    
    // -- parsing the CPU face file --
    void get_face_vector(  vector<double*>& face, vector<double> parsed );
    int parse_int( std::string line );
    double parse_double( std::string line, int* i, bool& success );
    vector<double> parse_vector( std::string line, int i );
    bool is_digit( char c );
    
    // -- if wall coupling is active --
    bool active;
    bool fileRead;
    
    
    int coord2bin(const double *x, int ic) const;

    // -- cpu's bounding box --
    double * coords_;  
    
    // -- whole domain bounding box --
    double * bounds__;
    
    int* n_send;
    int* n_recv;    
	
    // -- number of elements on CFD side --
    int n_cfd;
    
    // -- number of elements on the DEM side --
    int n_dem;
    
    int* localtagCFD_;
    
    // list CFD tasks where i:th element belongs
    std::vector<int> cpu_indexes;
    
    // -- element indices on DEM side being sent to CFD side --
    std::vector<int>* elementCPUMapping_;
    
    // -- element indices on CFD side being sent to DEM side --
    std::vector<int>* elementDEMMapping_;
    
    // -- local DEM indices of the CFD side array --
    std::vector<int> local_ids;
    
    // -- cpu ranks corresponding to CFD side array --
    std::vector<int> cpu_ids;  
    
    // -- routines for reconstructing elementDEMMapping
    
    virtual void send_data( int );	// int: destination
    virtual void recv_data( int );    // int: source
    virtual void self_comm();         // communication in own rank
    
    private :
    
    // -- DEM wall elements may lie slightly outside the bounding box computed from CFD faces, enlarge the bounding boces by small number --
    const static double asciiPrecission = 1e-4;
    
    // -- status object for nonblocking communication
    MPI_Request destReqs;
    MPI_Request send_reqs;
    MPI_Status send_stats;
   
};



}

#endif
