/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Efield
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nonspherical,FixNonSpherical)

#else

#ifndef LMP_FIX_NONSPHERICAL_H
#define LMP_FIX_NONSPHERICAL_H

#include "fix.h"
#include "math.h"

namespace LAMMPS_NS {
  
  class FixNonSpherical : public Fix 
  {
    
  public:
    FixNonSpherical(class LAMMPS *, int, char **);
    ~FixNonSpherical();
    
    virtual void post_create();
    virtual void pre_delete(bool unfixflag){ UNUSED(unfixflag); };
    
    virtual void setup_pre_exchange();
    
    virtual void pre_force(int vflag);
    virtual void post_force(int vflag);
    
    virtual void final_integrate();
    virtual int setmask();
    virtual void init();
    
  protected:    
    
    bool insertFlag;
    
    class PairGran *pair_gran;
    
    class FixPropertyParticle* fix_position;
    class FixPropertyParticle* fix_velocity;    
    class FixPropertyParticle* fix_normal;
    
    
    
    class FixPropertyParticle* fix_neighbor_index;    
    
    
    // -- list of particles that are neighbors --
    class FixPropertyParticle* fix_id;
    class FixPropertyParticle* fix_neighbor_id;
    
    double eqdistance;
    double kn;
    double kt;
    
    double gamman;
    double gammat;
    
    const static int neighborCount = 2;
    
    void initNeighbors();
    void integrateOrientation();
    
    bool isNeighbor(int,int) const;
    
    inline double getAngle( const double* a, const double* b ) const
    {
       double angle = acos( a[0] * b[0] + a[1] * b[1] + a[2] * b[2] );
    
       while( angle > 3.1415926/2.0 ) angle -= 3.1415926;
       while( angle < -3.1415926/2.0 ) angle += 3.1415926;
       
       return angle;
    }
    
    inline double dot( const double* a, const double* b ) const
    {
       return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
    
    inline void cross( const double* a, const double* b, double* c ) const
    {
	c[0] = a[1]*b[2] - b[1]*a[2];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
    }
    
    // -- check if particle is at the end of the particle chain --
    inline bool atEnd( int i ) const
    {
        
	for( int ii = 0; ii < neighborCount; ++ii )
	    if( fix_neighbor_id->array_atom[i][ii] < 0 ) 
	       return true;
	    	
	return false;
	
    }
    
    inline void appendToList( double* p, int i )
    {
        int index = 0;
	while( index < neighborCount && p[index] >= 0 )
	{
	   // -- i already in the list --
	   if( int( p[index] ) == i ) return;
	   
	   ++index;
	}
	if( index < neighborCount ) p[index] = double(i);
	
    }
    
    inline int listLength( double* p ) const
    {
        
	for( int i = 0; i < neighborCount; ++i )
	    if( p[i] < 0 ) return i;
	
	return neighborCount;
	
    }
    
    inline double mag( const double* a ) const
    {
        return sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
    }
    
    inline void normalize( double* a )
    {
       const double magA = sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
       for( int i = 0; i < 3; ++i )
          a[i] /= magA;
    }
    
    inline double abs( double a )
    {
       return a > 0 ? a : -a;
    }
            
    
  };  

}

#endif
#endif
