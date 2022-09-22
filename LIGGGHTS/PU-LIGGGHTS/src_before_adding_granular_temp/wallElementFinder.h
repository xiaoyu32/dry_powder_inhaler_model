/*
 * Jari Kolehmainen, 2019 
 * Determines if a given point lies on the CFD bounding mesh
 */


#ifndef _WALL_ELEMENT_FINDER_H_
#define _WALL_ELEMENT_FINDER_H_

#include <vector>
#include <cstdio>

class WallElementFinder
{

  protected :

   class Face
   {
   
      public :
       
       Face( char* byteStream );
       
       Face( const Face& f );
       Face( const std::vector<double*>& f );      
       Face();
       
       ~Face();
       
       void getBounds( double*, double* ) const;
       
       void addPoint( double* v );
       
       bool isOnFace( const double* v ) const;
       
       void extractByteStream( char* ) const;
       
       int byteSize() const;
       
       inline int npoints() const
       {
           return points.size();
       }
       
       void printPoints() const;
       
      protected :
       
       // -- precomputes orthogonal basis matrix for the plane spanned by the face, and projects all the faces on the plane --
       void init();	
	
       // -- compute counter clockwise angle of vectors a and b --
       double vectorAngle( const double* a, const double* b ) const;
	
       double mag( const double* a ) const;	
       double dot( const double* a, const double* b ) const;	
	
	// -- inline functions --
       inline double abs( const double& a ) const
       {
	   return a > 0 ? a : -a;
       }
       
       inline void normalize( double* a ) const
       {
           const double magA = mag( a );
	   for( int i = 0; i < 3; ++i )
	      a[i] /= magA;
       }
       
       // -- compute (b-a) & e --
       double projectDistance( const double* a, const double* b, const double* e ) const;  
       
       void projectPoint( const double* a, const double* b, double* c ) const;
       
       inline void minus( const double* a, const double* b, double* c ) const
       {
	   for( int i = 0; i < 3; ++i )
	      c[i] = b[i]-a[i];       
       }
       
       inline void cross( const double* a, const double* b, double* c ) const
       {
	   c[0] = a[1]*b[2] - a[2]*b[1];
	   c[1] = a[2]*b[0] - a[0]*b[2];
	   c[2] = a[0]*b[1] - a[1]*b[0];          
       }
       
       inline void cpyVector( const double* a, double* b ) const
       {
           for( int i = 0; i < 3; ++i )
	      b[i] = a[i];
       }
       
       bool initFlag;
       bool initCalled;
       
       double* refPoint;
       double basisAngle;
       
       
       // -- orthogonal basis vector; first two vectors span the plane and third one is orthogonal to the face --
       double Q[3][3];
       
       // -- face edge points projected on the face plane by Q --
       std::vector<double*> projectedPoints;
       
       // -- face points --
       std::vector<double*> points;
       
       // -- maximum face deviance from plane --
       double distExtra;
       
       // -- maximum allowed distance from the plane spanned by the faces (should be ASCII %f accuracy) --
       static const double distanceThreshold = 1.0e-6;

       // -- maximum allowed difference from the full angle (in radians) --
       static const double angleThreshold = 1e-2;

       // -- full angle in radians --
       static const double fullAngle = 2.0 * 3.1415926;       
	
       // -- large number --
       static const double INF = 1e9;	
       
       // -- minimum angle between face edge vectors before take can be used to span base --
       static const double minAngle = 0.25;
	
   };

   // -- return index of the octoTree member, which contains face -- 
   int containedIndex( const WallElementFinder::Face& face ) const;

   // -- return index of the sub-container, -1 if point is not in the subcontainer --
   int pointIndex( const double* point ) const;
   
   bool addNewFace( const WallElementFinder::Face& ff );
   
   void set( double*, double* );
   
   void byteStream( char*, int& ) const;
   
   // -- container bounds --
   double lower[3];
   double upper[3];
   double mid[3];   
   
   // -- number of elements --
   int nelements;
   
   std::vector<WallElementFinder::Face> faces; 
   WallElementFinder* octoTree[8]; 
   
   bool isContained( const WallElementFinder::Face& face ) const;
       
  public :
    
    // -- construct from bounds --
    WallElementFinder(double*, double*);
    
    // -- construct from byte stream --
    void initFinder( char* stream, int bytes );
    
    void clearAll();
    
    ~WallElementFinder();

    bool addFace( std::vector<double*>& face_ );
    
    void printBounds();
    
    bool isOnFaces( const double* ) const;
    
    inline int numberOfElements() const
    {
       return nelements;
    } 
    
    inline int nfaces() const
    {
       return faces.size();
    }
    
    char* extractByteStream() const;
    
    int byteSize() const;
    
};




#endif
