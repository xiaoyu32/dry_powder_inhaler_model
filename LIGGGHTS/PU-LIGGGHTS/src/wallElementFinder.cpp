/*
 * Jari Kolehmainen, 2019 
 * Determines if a given point lies on the CFD bounding mesh
 */

#include "wallElementFinder.h"
#include "math.h"

// -- face class members --

WallElementFinder::Face::Face() :
initFlag( false ),
initCalled( false ),
distExtra( 0 ),
refPoint( NULL ),
basisAngle( 0 )
{}

WallElementFinder::Face::Face( const WallElementFinder::Face& f ):
initFlag( false ),
initCalled( false ),
distExtra( 0 ),
basisAngle( 0 )
{
   for( int i = 0; i < f.points.size(); ++i )
       this->addPoint( f.points.at(i) );
}

WallElementFinder::Face::Face( char* byteStream ) :
initFlag( false ),
initCalled( false ),
distExtra( 0 ),
refPoint( NULL ),
basisAngle( 0 )
{
   
    int* nPtr = reinterpret_cast<int*>( byteStream );
    int nn = nPtr[0]; 
    
    double* ptr = reinterpret_cast<double*>( byteStream + sizeof( int ) );
          
    for( int i = 0; i < nn; ++i )
    {
	this->addPoint( &ptr[3*i] );
    }

     
}

WallElementFinder::Face::Face( const std::vector<double*>& f ):
initFlag( false ),
initCalled( false ),
distExtra( 0 ),
refPoint( NULL ),
basisAngle( 0 )
{
   for( int i = 0; i < f.size(); ++i )
       this->addPoint( f.at(i) );    
}

WallElementFinder::Face::~Face()
{
   
   for( int i = 0; i < points.size(); ++i )
      delete[] points.at(i);
   
   points.clear();
   
   if( initCalled )
   {
      for( int i = 0; i < projectedPoints.size(); ++i )
          delete[] projectedPoints.at(i);
      projectedPoints.clear();
   }
   
}


int WallElementFinder::Face::byteSize() const
{
    int byteSize_ = 0;
    byteSize_ += sizeof( int ); // -- number of points in face --
    
    // -- face points --
    for( int i = 0; i < points.size(); ++i )
       byteSize_ += 3*sizeof( double );
    
    return byteSize_;
  
}

void WallElementFinder::Face::extractByteStream( char* stream ) const
{
        
    int* nPtr = reinterpret_cast<int*>(stream); 
    nPtr[0] = points.size();
    
    double* ptr = reinterpret_cast<double*>( stream + sizeof( int ) );
    
    int index = 0;
    for( int i = 0; i < points.size(); ++i )
    {
        double* v = points.at(i);
        for( int j = 0; j < 3; ++j )
	    ptr[index++] = v[j];
	
    }
    
}

void WallElementFinder::Face::init()
{
    
    if( !initFlag ) return; // -- too few points to initialize --
    
    bool basisChanged = false;
    
    if( !initCalled || basisAngle < minAngle ) // -- initialize the otrhogonal base matrix --
    {

       // -- spanning vectors of the plane where face belongs --
       double e1[3] = {0,0,0};
       double e2[3] = {0,0,0};
       double e3[3] = {0,0,0};
           
       int index = 2;	   
        
	      
       do
       {
       	  int modulo = points.size();

	  // -- face must have at least three points --
	  double* p1 = points.at( 0 );
	  double* p2 = points.at( 1 );
	  double* p3 = points.at( index );   
	  
	  refPoint = p1;
	  
	  minus( p1,p2,e1 );
	  minus( p1,p3,e2 );

	  normalize( e1 );
	  normalize( e2 );
	  cross( e1, e2, e3 );
	  
	  this->basisAngle = this->abs( vectorAngle( e1,e2 ) );
	  ++index;
	  
       }while( 
		this->basisAngle <= minAngle
       		&& index < points.size() 
	     );
       
       basisAngle = this->abs( vectorAngle( e1,e2 ) );
       
       normalize( e3 );

       // -- create orthogonal basis for the plane --
       double ee2[3] = {0,0,0};
       cross( e1,e3,ee2 );

       this->Q[0][0] = e1[0];
       this->Q[0][1] = e1[1];
       this->Q[0][2] = e1[2];

       this->Q[1][0] = ee2[0];
       this->Q[1][1] = ee2[1];
       this->Q[1][2] = ee2[2];

       this->Q[2][0] = e3[0];
       this->Q[2][1] = e3[1];
       this->Q[2][2] = e3[2];    

       initCalled = true;
       basisChanged = true;       
    }

    double distExtra_ = 0;
    
    if( basisChanged ) // -- basis has changed recalculate the projected points --
    {
       this-> distExtra = 0;
       
       // --check if angle threshold is reached --
       for( int i = 0; i < projectedPoints.size(); ++i )
       {
 	   delete[] projectedPoints.at(i);
       }    
       
       projectedPoints.clear();
       
    }
    
    // -- use the new points to determine the maximum deviance of the face from plane --
    for( int i = projectedPoints.size(); i < points.size(); ++ i )
    {
        distExtra_ = this->abs( projectDistance( refPoint, points.at(i), this->Q[2] ) );	
	if( distExtra_ > distExtra ) distExtra = distExtra_;
    }
      
    // -- project new points --    
    for( int i = projectedPoints.size(); i < points.size(); ++ i )
    {
        double* cc = new double[3];
    	projectPoint( refPoint, points.at(i), cc );
	
	// -- remove numerical remnants from the projected points (should lie in the plane) --
	projectedPoints.push_back( cc );
	
    }

    
}

void WallElementFinder::Face::printPoints() const
{
    
    for( int ii = 0; ii < points.size(); ++ii )
    {
        printf( "%d : [%e %e %e] \n", ii, points.at(ii)[0], points.at(ii)[1], points.at(ii)[2] );
    }
    
}

bool WallElementFinder::Face::isOnFace( const double* v ) const
{
    
    //printf( "isOnFace Called ! \n" );
    
    // -- face is incomplete --
    if( !initCalled ) return false;

    if( basisAngle < minAngle )
    {
       printf( "Warning: small base vector angle detected. Basis might be ill-conditioned!\n" );
    }

    
    const double* e3 = this->Q[2];
    
    // -- distance from the place --
    const double dist = projectDistance(  points.at(0), v, e3 );
    
    // -- point is not in the place --
    if( this->abs( dist ) > distanceThreshold + this->distExtra ) return false;


    // -- check if the point is on the face --


    // -- map the face and the point of interest on the plane --
    double projected_v[3] = {0,0,0};
    projectPoint( refPoint, v, projected_v );
    
    //printf( "projected_v = [ %e %e %e ] \n", projected_v[0], projected_v[1], projected_v[2] );
    
    double angleSum = 0;

    // -- start traversing the face from the last point --
    double pp[3] = {0,0,0};
    double pn[3] = {0,0,0};
    
    minus( projected_v, projectedPoints.at( projectedPoints.size() - 1 ), pp );
    
    //printf( "projectedPoints.size() = %d   points.size() = %d \n", projectedPoints.size(), points.size() );
    
    for( int i = 0; i < projectedPoints.size(); ++i )
    {
    
        minus( projected_v, projectedPoints.at(i), pn );
	
	//printf( "projectedPoints.at(i) = [ %e %e %e ] \n", projectedPoints.at(i)[0], projectedPoints.at(i)[1], projectedPoints.at(i)[2] );
	
	const double vangle = vectorAngle( pp, pn );
		
	//printf( "Angle %d  =  %f \n ", i, 360.0 * vangle/(fullAngle) );
	
	// -- add angle --
	angleSum += vangle;
	cpyVector(pn,pp);
	
    }
        
    //printf( "angleSum = %e \n", 360.0 * angleSum/fullAngle );
    
    // -- if angles sum up to full round point point is on the face --
    if( 
    	this->abs( angleSum - fullAngle ) < angleThreshold || /* -- counter clockwise -- */
	this->abs( angleSum + fullAngle ) < angleThreshold    /* -- clockwise -- */
      )
    {
	return true;
    }else
    {
	return false;
    }

}

void WallElementFinder::Face::addPoint( double* v )
{
    
    double* v_ = new double[3];
    for( int i = 0; i < 3; ++i )
       v_[i] = v[i];
       
    points.push_back( v_ );
    
    if( points.size() >= 3 ) // -- there are at last three points in face -> initialize face projection --
    {
       this->initFlag = true;
       this->init();
    }
    
}

void WallElementFinder::Face::getBounds( double* min, double* max ) const
{

    min[0] = INF;
    min[1] = INF;
    min[2] = INF;

    max[0] = -INF;
    max[1] = -INF;
    max[2] = -INF;	   

    for( int i = 0; i < points.size(); ++i )
    {
	double* p = points.at(i);

	for( int ii = 0; ii < 3; ++ii )
	{
	    if( p[ii] < min[ii] ) min[ii] = p[ii];
	    if( p[ii] > max[ii] ) max[ii] = p[ii];
	}
    }
    
}

// compute (b-a) & e
double WallElementFinder::Face::projectDistance( const double* a, const double* b, const double* e ) const
{
   double c[3];
   minus( a,b,c );
   return dot( c, e );
}

void WallElementFinder::Face::projectPoint( const double* a, const double* b, double* c ) const
{

    for( int i = 0; i < 3; ++i )
    {   
       c[i] = 0;
       for( int j = 0; j < 3; ++j )
       {
           c[i] += this->Q[i][j] * (b[j] - a[j]);
       }
    }   
    
    // -- set third component zero --
    c[2] = 0.0;   
}

double WallElementFinder::Face::vectorAngle( const double* a, const double* b ) const
{
    double out[3] = {0,0,0}; 
    
    cross( a,b,out );
    const double in = dot(a,b)/( mag(a) * mag(b) );
        
    if( out[2] < 0 ) // -- counter clockwise --
    {
        return acos( in );
    }else{ // -- clockwise --
        return - acos( in );
    }
}

double WallElementFinder::Face::mag( const double* a ) const
{
   
   double res = 0;
   for( int i = 0; i < 3; ++i )
      res += a[i]*a[i];
   return sqrt( res );
   
}

double WallElementFinder::Face::dot( const double* a, const double* b ) const
{
   double res = 0;
   for( int i = 0; i < 3; ++i )
      res += a[i]*b[i];
   return res;
}


// -- wall element finder class --

WallElementFinder::WallElementFinder( double* lower_, double* upper_ ) :
nelements( 0 )
{
   set( lower_, upper_ );
   for( int i = 0; i < 8; ++i )
       octoTree[i] = NULL;
}

WallElementFinder::~WallElementFinder()
{
   faces.clear();
   for( int i = 8; i < 8; ++i )
      if( octoTree[i] ) delete octoTree[i];
}

void WallElementFinder::clearAll()
{
   this->nelements = 0;
   faces.clear();
   for( int i = 8; i < 8; ++i )
      if( octoTree[i] ) octoTree[i]->clearAll();   
}


int WallElementFinder::containedIndex( const WallElementFinder::Face& face ) const
{

    double min[3];
    double max[3];

    face.getBounds( min, max );

    int ii[3] = {-1,-1,-1};

    for( int i = 0; i < 3; ++i )
    {

	if( max[i] < this->mid[i] )
	{   
	   ii[i] = 0;
	   continue;
	}

	if( min[i] >= this->mid[i] )
	{  
	   ii[i] = 1;
	   continue;  
	}   

	return -1;

    }

    return ii[0] + ii[1] * 2 + ii[2] * 4;

}

int WallElementFinder::pointIndex( const double* point ) const
{
   int ii[3] = {-1,-1,-1};

   for( int i = 0; i < 3; ++i )
   {
      if( point[i] < lower[i] || point[i] > upper[i] ) return -1;

      if( point[i] < mid[i] )
      {  
	  ii[i] = 0;
	  continue;
      }

      if( point[i] >= mid[i] )
      {  
	  ii[i] = 1;
	  continue;
      }
   }

   return ii[0] + ii[1] * 2 + ii[2] * 4;
}

// -- check if the face is fully enclosed in the container --
bool WallElementFinder::isContained( const WallElementFinder::Face& face ) const
{

    double min[3];
    double max[3];

    face.getBounds( min, max );

    for( int ii = 0 ; ii < 3; ++ii )
    {

	if( min[ii] > lower[ii] && max[ii] < upper[ii] ) continue;

	return false;

    }

    return true;

}

bool WallElementFinder::addNewFace( const WallElementFinder::Face& ff )
{
    
    // -- check if face is contained in any of the sub-containers --
    int index = containedIndex( ff );
        
    if( index >= 0 )
    {

	if( octoTree[index] ) return octoTree[index]->addNewFace( ff );
	else
	{
	   // -- create new sub-container --
	   int ii[3] = {0,0,0};
	   int index_ = index;


	   ii[0] = index_%2;

	   index_ = ( index_ - ii[0] )/2;
	   ii[1] = index_%2;
	   index_ = ( index_ - ii[1] )/2;

	   ii[2] = index_%2;

	   double newLower[3];
	   double newUpper[3];

	   for( int i = 0; i < 3; ++i )
	   {
	       newLower[i] = ii[i] == 0 ? lower[i] : mid[i];
	       newUpper[i] = ii[i] == 0 ? mid[i] : upper[i];
	   }

	   octoTree[index] = new WallElementFinder( newLower, newUpper );
	   return octoTree[index]->addNewFace( ff );

	}

    }else // -- face is not contained in the sub-containers -> add face in the present container
    {
	faces.push_back( ff );
	return true;
    }
        
    // -- should never come here --
    return false;

}

void WallElementFinder::set( double* lower_, double* upper_ )
{
    for( int ii = 0; ii < 3; ++ii )
    {
       lower[ii] = lower_[ii];
       upper[ii] = upper_[ii];
       mid[ii] = (upper[ii] + lower[ii])/2.0; 
    }
}

bool WallElementFinder::isOnFaces( const double* point ) const
{

    int index = pointIndex( point );

    // -- check if point is on the faces in the container --

    for( int i = 0; i < faces.size(); ++i )
    {
	const WallElementFinder::Face& ff = faces.at(i);

	if( ff.isOnFace( point ) ) // -- face found -- 
	{		     
	    return true;
	}
    }

    if( index < 0 ) return false;

    if( octoTree[index] )
    {
	return octoTree[index]->isOnFaces( point );
    }	     
    else return false;
    
}

char* WallElementFinder::extractByteStream() const
{
    int index = 0;
    
    char* stream = new char[ this->byteSize() ];
    this->byteStream( stream, index );
    
    return stream;
    
}

void WallElementFinder::byteStream( char* stream, int& index ) const
{
    for( int i = 0; i < faces.size(); ++i )
    {
        faces.at(i).extractByteStream( stream + index );
	index += faces.at(i).byteSize();
    }
    
    for( int i = 0; i < 8; ++i )
    {
        if( octoTree[i] )
	{
	    octoTree[i]->byteStream(stream,index);
	}
    }    
    
}

int WallElementFinder::byteSize() const
{
        
    int byteSize_ = 0;
    
    for( int i = 0; i < faces.size(); ++i ) // -- size of faces belonging to this finder --
    {
        byteSize_ += faces.at(i).byteSize();
    }
    
    for( int i = 0; i < 8; ++i )
    {
        if( octoTree[i] )
	{
	    byteSize_ += octoTree[i]->byteSize();
	}
    }

    return byteSize_;
    
}

void WallElementFinder::initFinder( char* stream, int bytes )
{
    
    int index = 0;
        
    while( index < bytes )
    {
        Face ff( stream + index );	
        index += ff.byteSize();
	
	if( isContained( ff ) )
	{
	    ++nelements;   
	    this->addNewFace( ff );
	}else{
	    printf( "Failed to add face \n" );
	    ff.printPoints();
	}
    }
    
}

bool WallElementFinder::addFace( std::vector<double*>& face_ )
{

   WallElementFinder::Face ff( face_ );

   if( !isContained( ff ) )
   {
      /*
      printf( " \n" );
      printf( "lower = [ %e %e %e ] \n", lower[0], lower[1], lower[2] );
      printf( "upper = [ %e %e %e ] \n", upper[0], upper[1], upper[2] );
      ff.printPoints();
      printf( " \n" );
      */

      return false;
   }

   ++nelements;       
   return this->addNewFace( ff );

}

void WallElementFinder::printBounds()
{
   printf( "lower = [ %e %e %e ] \n", lower[0], lower[1], lower[2] );
   printf( "upper = [ %e %e %e ] \n", upper[0], upper[1], upper[2] );
}





