#ifndef _KD_TREE_SEARCH_
#define _KD_TREE_SEARCH_

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <map>
#include "math.h"

class KdTreeSearch
{
    
    public :
    
    static const int SMALLER = -1;
    static const int EVEN = 0;
    static const int LARGER = 1;
    static const int ENCLOSE = 2;
    
    static const double _small_ = 1e-13;
    
    double* dir[2];
    
    class Node;
    
    // -- list of nodes in the tree --
    std::vector<Node*> nodeList; 
    
    // -- map containing ids of Nodes that have already been intersected (or cannot be intersected) by a ray --
    mutable std::map<int,bool> intersected;
    int id_counter;
    
    
    // -- root of the search tree --
    Node* root;
    
    KdTreeSearch( double* dir1_, double* dir2_ ) 
    {
        dir[0] = dir1_;
	dir[1] = dir2_;
	
        root = NULL;
	id_counter = 0;
    }
    
    ~KdTreeSearch()
    {
       if( root ) delete root;
    }
    
    
    class Node
    {
        
	public :
	
	const int id;
	
	Node* lower;
	Node* upper;
	Node* midSmall;
	Node* midLarge;
	
	Node* parent;
	
	// -- list of neighboring nodes --
	std::vector<Node*> neighbors;
		
	const std::vector< double* >& face;
	
	Node( const std::vector< double* >& face_, int id_ ) :
	face( face_ ),
	id( id_ )
	{
	   
	   lower = NULL;
	   upper = NULL;
	   midSmall = NULL;
	   midLarge = NULL;
	   
	   parent = NULL;
	   
	};
	
	~Node()
	{
	   if( lower ) delete lower;
	   if( midSmall ) delete midSmall;
	   if( midLarge ) delete midLarge;
	   if( upper ) delete upper;
	}
	
	void addNeighbor( Node* n )
	{
	    
	    if( this == n ) return;
	    
	    for( int i = 0; i < neighbors.size(); ++i )
	        if( n == neighbors[i] ) return;
	    
	    for( int i = 0; i < face.size(); ++i )
	    {
	        for( int j = 0; j < n->face.size(); ++j )
		{
		    if( distance( face.at(i), n->face.at(j) ) < KdTreeSearch::_small_ )
		    {
		        neighbors.push_back( n );
			n->neighbors.push_back( this );
			return;
		    }
		}
		
	    }
	    
	} 
	
	void printFace() const
	{
	    std::cout<<"Face ("<<id<<") ";
	    for( int i = 0; i < face.size(); ++i )
	       std::cout<<"( "<<face.at(i)[0]<<", "<<face.at(i)[1]<<", "<<face.at(i)[2]<<" ) ";
	    std::cout<<std::endl;
	}

	
	inline double distance( const double* a, const double* b ) const
	{
	    double rr = 0;
	    
	    for( int i = 0; i < 3; ++i )
	       rr += ( a[i]-b[i] ) * ( a[i]-b[i] );
	    return sqrt( rr );   
	}
	
	double dotProduct( const double* a, const double* b ) const
	{
	     double r = 0;
	     for( int ii = 0; ii < 3; ++ii )
	        r += a[ii] * b[ii];
	     
	     return r;
	     
	}
	
	void getMinMax( const std::vector< double* >& ff, const double* dir, double& min, double& max ) const
	{
	    
	    double d = dotProduct( ff[0], dir );
	    
	    min = d;
	    max = d;
	    
	    for( int i = 1; i < ff.size(); ++i )
	    {
	        // coordinate in the dir direction
		d = dotProduct( ff[i], dir );
		
		if( d < min ) min = d;
		else if( d > max ) max = d;
		
		
	    }
	    
	}
	
	int compare( const std::vector< double* >& newFace, const double* dir ) const
	{
	    
	    double min1;
	    double max1;
	    
	    double min2;
	    double max2;	    
	    
	    getMinMax( face, dir, min1, max1 );
	    getMinMax( newFace, dir, min2, max2 );

	    if( min2 < min1 && max2 < max1 ) return SMALLER;
	    if( max2 > max1 && min2 > min1 ) return LARGER;
	    
	    if( min2 <= min1 && max2 >= max1 ) return ENCLOSE;
	    
	    return EVEN;
	    
	}
	
	int compare( const double* point, const double* dir ) const
	{
	    
	    double min;
	    double max;	    
	    getMinMax( face, dir, min, max );
	    
	    double projection = dotProduct( point, dir );
	    
	    if( projection < min ) return SMALLER;
	    if( projection > max ) return LARGER;
	    
	    return EVEN;
	    
	}
	
	
    }; // -- end of class node --
    
    void addFace( const std::vector< double* >& face )
    {

	if( !root )
	{ 
	    root = new Node( face, id_counter++ );
	    nodeList.push_back( root );
	}
	else
	{ 
	    addFaceToNode( face, root, 0 );   
	}
	
    }

    void initNeighborList()
    {
        // fixme: n^2 loop, do something fast (done only once per simulation, so maybe acceptable)
	for( int i = 0; i < nodeList.size(); ++i )
	{  
	    Node* n = nodeList.at(i); 
	     
            for( int j = i+1; j < nodeList.size(); ++j )
	    {
	        n->addNeighbor( nodeList.at(j) );
	    }
	}
	
    }

    void initTree()
    {
	
	initNeighborList();
	
	/*std::cout<<"Number of neighbors: "<<nodeList.size()<<std::endl;
	
	for( int i = 0; i < nodeList.size(); ++i )
	{
	    std::cout<<"Node: "<<i<<" neighbors: "<<nodeList.at(i)->neighbors.size()<<std::endl;
	}*/
	
    }
     
    
    void addFaceToNode( const std::vector< double* >& face, KdTreeSearch::Node* n, int dirCounter )
    {
        
	int cc = n->compare( face, dir[dirCounter%2] );
		
	if( cc == SMALLER )
	{
	    
	    if( !n->lower )
	    {
	        n->lower = new Node( face, id_counter++ );
		n->lower->parent = n;
		
		nodeList.push_back( n->lower );
	    }else
	    {
	        addFaceToNode( face, n->lower, (dirCounter+1)%2 );
	    }
	    
	}else if( cc == LARGER )
	{

	    if( !n->upper )
	    {
	        n->upper = new Node( face, id_counter++ );
		n->upper->parent = n;
		
		nodeList.push_back( n->upper );
	    }else
	    {
	        addFaceToNode( face, n->upper, (dirCounter+1)%2 );
	    }
	
	}else if( cc == EVEN ) // new face is fully enclosed in the present face
	{
	    if( !n->midSmall )
	    {
	        n->midSmall = new Node( face, id_counter++ );
		n->midSmall->parent = n;
		
		nodeList.push_back( n->midSmall );
	    }else
	    {
	        addFaceToNode( face, n->midSmall, (dirCounter+1)%2 );
	    }
	}else if( cc == ENCLOSE ) // new face fully encloses the present face i.e. is larger
	{
	    if( !n->midLarge )
	    {
	        n->midLarge = new Node( face, id_counter++ );
		n->midLarge->parent = n;
		
		nodeList.push_back( n->midLarge );
	    }else
	    {
	        addFaceToNode( face, n->midLarge, (dirCounter+1)%2 );
	    }
	}
	
    }    
    
    double abs( const double a ) const
    {
        return a>0 ? a : -a;
    }
    
    int intersectionPoints( const double* RESTRICT ray, const double* RESTRICT point ) const
    {
        // -- clear intersected points --
	intersected.clear();
	
	//std::cout<<"BEGIN"<<std::endl;
	//std::cout<<std::endl;
	
	if( !root ){
	   return 0;
	}else
	{
	
	    int numberIntersections = 0;
	    intersectionPoint( numberIntersections, root, ray, point, 0 );
	    
	    /*
	    std::cout<<"intersectionPoints: "<<numberIntersections<<std::endl;
	    
	    for( std::map<int,bool>::iterator it = intersected.begin(); it!=intersected.end(); ++it)
	    {
	        std::cout<<it->first<<" ";
	        nodeList.at( it->first )->printFace(); 
	    }
	    std::cout<<std::endl<<std::endl;
	    */
	    
	    return numberIntersections;
	    
	}
	
    }
    
    void normalize( double* c ) const
    {
        
	double r = 0;
	for( int i = 0; i < 3; ++i )
	    r += c[i]*c[i];
	
	r = sqrt(r);
	
	for( int i = 0; i < 3; ++i )
	   c[i] /= r;

	
    }
    
    void crossProduct( const double* RESTRICT a, const double* RESTRICT b, double* RESTRICT c ) const
    {
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];        
    }

    double dotProduct( const double* a, const double* b ) const
    {
	 double r = 0;
	 for( int ii = 0; ii < 3; ++ii )
	    r += a[ii] * b[ii];

	 return r;

    }
    
    void findBounds( 	const double* p, const double* a, const double* b, const double* l, 
    			double& min, double& max ) const
    {
        
	double t[3];
	
	t[0] = dotProduct( p, l );
	
	double p2[3];
	double p3[3];
	
	for( int i = 0; i < 3; ++i )
	{
	    p2[i] = p[i] + a[i];
	    p3[i] = p[i] + b[i];
	}
	
	t[1] = dotProduct( p2, l );
	t[2] = dotProduct( p3, l );
	
	min = t[0];
	max = t[0];
	
	
	for( int i = 1; i < 3; ++i )
	{
	    
	    if( t[i] < min ) min = t[i];
	    if( t[i] > max ) max = t[i];
	
	}

    }
    
    // return -1 if not intersected, 1 if intersected in the interior, 2 if intersected on one of the boundaries
    int solveTriangle( const double* RESTRICT pt, const double* RESTRICT a, const double* RESTRICT b, const double* RESTRICT p, const double* RESTRICT ray ) const
    {
	
	
	
	//compute determinate of matrix [a b -ray]   
	double det = determinant( a, b, ray );

	// - ray is parallel to the triangle plane -
	if( abs( det ) < 1e-15 )
	{

	    //DO SOMETHING HERE (voting should take care of this even when point lies in the plane)
	    //std::cout<<"Ray = "<<ray[0]<<" "<<ray[1]<<" "<<ray[2]<<std::endl;
	    //std::cout<<"Zero determinant..."<<std::endl;
	    /*
	    double c[3];
	    
	    // - get orthogonal direction to the triangle element - 
	    crossProduct(a,b,c);
	    normalize( c );
	    
	    double d1 = dotProduct( c, p );
	    double d2 = dotProduct( c, pt );
	    
	    if( abs( d1-d2 ) < 1e-9 ) //particle center lies in the same plane as the triangle
	    {
	        
		//check if the ray intersect triangle by looking at the interval projected in the orthogonal direction of the ray
	        double cc[3];
		
		crossProduct(ray, c, cc );
	    	normalize( cc );
		
		double q1 = dotProduct( cc, pt );
		
		double min;
		double max;
		
		findBounds( p, a, b, cc, min, max );
		
		if( q1 > min && q1 < max ) return true;
		else			   return false;
		
	    }else{ // points do not lie in the same plane spanned by triangle and hence cannot intersect
	       return false;
	    }
	    */
	    
	    return -1;
	}


	//there will be intersection with the triangle plane and a line going through the point
	double sol[3];

	double rhs[3];

	for( int i = 0; i < 3; ++i )
	{
            rhs[i] = p[i] - pt[i];
	}

	gaussianSolve( sol, a, b, ray, rhs );

        // -- intersecting at the edge of triangle --
	if( 
	    ( sol[0] >= -KdTreeSearch::_small_ && sol[0] <= 0 ) || 
	    ( sol[0] >= -KdTreeSearch::_small_ && sol[0] <= 0 ) || 
	    ( sol[0] + sol[1] >= 1.0 && sol[0] + sol[1] < 1.0 + KdTreeSearch::_small_ ) 
	  )
	{
	    // -- handle seprately --
	    return 2;
	}

	//if the intersection point lies within the triangle
	if( sol[2] < 0 && sol[0] >= 0 && sol[1] >= 0  && (sol[0] + sol[1] <= 1.0 ) )
	{
	    //std::cout<<"determinant = "<<det<<std::endl;
	    //std::cout<<"sol = "<<sol[0]<<" "<<sol[1]<<" "<<sol[2]<<std::endl;
	    //std::cout<<"p = "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
	    //std::cout<<"ray = "<<ray[0]<<" "<<ray[1]<<" "<<ray[2]<<std::endl;
	    //std::cout<<std::endl;
            return 1;
	}else{
            return -1;
	}

    } 

    double determinant( const double* a, const double* b , const double* c ) const
    {

	return   a[0] * ( b[1] * c[2] - c[1] * b[2] ) 
    	       - a[1] * ( b[0] * c[2] - c[0] * b[2] )
    	       + a[2] * ( b[0] * c[1] - c[0] * b[1] );
    }


    void swap( double* a, double* b, int i ) const
    {

	if( abs( a[i] ) < abs( b[i] ) )
	{

	    double buf[4] = {b[0], b[1], b[2], b[3]};

	    for( int k = i; k < 4; ++k )
	    {
		b[k] = a[k];
		a[k] = buf[k];
	    }

	} 

    }

    void sum( const double* a, double* b, int i ) const
    {

	if( abs( a[i] ) < 1e-16 ){

	    //this should never happen
	    std::cout<<"Gaussian Elimination Failed !!!"<<std::endl;

            return;
	}

	const double cof = -b[i]/a[i];

	for( int k = i; k < 4; ++k )
	{

	    b[k] += cof* a[k];

	}

    } 

    void gaussianSolve( double* RESTRICT sol, const double* RESTRICT a , const double* RESTRICT b, const double* RESTRICT c, const double* RESTRICT rhs ) const
    { 

	double A[3][4] = {
    			    0,0,0,0,
			    0,0,0,0,
			    0,0,0,0
			 };

	//cpy to matrix
	for( int k = 0; k < 3; ++k )
	{   
	   A[k][0] = a[k];
	   A[k][1] = b[k];
	   A[k][2] = c[k];

	   A[k][3] = rhs[k];
	}

	//LU decompositions
	for( int i = 0; i < 3; ++i )
	{
	    //esure that the i:th leading element is non-zero
	    for( int j = i+1; j <3; ++j )
		swap( A[i], A[j], i );        

	    //add i:th line to other lines to zero i:th row components below i:th line
	    for( int j = i+1; j < 3; ++j )
		sum( A[i], A[j], i );

	}

	//back substitution
	for( int i = 2; i >= 0; --i )
	{

	    double res = A[i][3];

	    for( int j = i+1; j < 3; ++j )
	    {
		res -= A[i][j] * sol[j];
	    }

	    sol[i] = res/A[i][i];

	}

    }
    
    void intersectionPoint( int& numberIntersections, const Node* n, const double* RESTRICT ray, const double* RESTRICT point, int dirCounter ) const
    {
        
	int cc = n->compare( point, dir[dirCounter%2] );
		
	if( cc == SMALLER )
	{
	    
	    if( n->lower )
	    {
	        intersectionPoint( numberIntersections, n->lower, ray, point, (dirCounter+1)%2 );
	    }

	    if( n->midLarge )
	    {
	        intersectionPoint( numberIntersections, n->midLarge, ray, point, (dirCounter+1)%2 );
	    }
	    
	}else if( cc == LARGER )
	{

	    if( n->upper )
	    {
	        intersectionPoint( numberIntersections, n->upper, ray, point, (dirCounter+1)%2 );
	    }

	    if( n->midLarge )
	    {
	        intersectionPoint( numberIntersections, n->midLarge, ray, point, (dirCounter+1)%2 );
	    }
	
	}else if( cc == EVEN )
	{
	    // -- check face if it has at least 3 corners and its neighbors are not intersected --
	    if( n->face.size() >= 3 && intersected.count( n->id ) == 0 )
	    {
	        		
		const double* p1 = n->face[0];
		double p2[3];
		double p3[3];

		for( int i = 2; i < n->face.size(); ++i )
		{

		    for( int k = 0; k < 3; ++k )
		    {
			p2[k] = n->face[i-1][k] - p1[k];
			p3[k] = n->face[i][k] - p1[k];
		    }
		    
		    // ray can intersect single face atmost once
		    int triangleRes = solveTriangle( p1, p2, p3, point, ray );
		    
		    
		    if( triangleRes > 0 )
		    {
		       //std::cout<<"P1 = "<<n->face[0][0]<<" "<<n->face[0][1]<<" "<<n->face[0][2]<<std::endl;
		       //std::cout<<"P2 = "<<n->face[i-1][0]<<" "<<n->face[i-1][1]<<" "<<n->face[i-1][2]<<std::endl;
		       //std::cout<<"P3 = "<<n->face[i][0]<<" "<<n->face[i][1]<<" "<<n->face[i][2]<<std::endl;
		       
		      // n->printFace();
		       //std::cout<<std::endl;
		       // -- mark all the neighbors intersected --
		       if( triangleRes == 2 )
		       {
		           //std::cout<<"triangleRes == 2 "<<std::endl;
			   for( int j = 0; j < n->neighbors.size(); ++j )
			   {
		               intersected.insert( std::pair<int,bool>( n->neighbors.at(i)->id, true ) );
			   }
		       }
		       
		       intersected.insert( std::pair<int,bool>( n->id, true ) );
		       
		       ++numberIntersections;
		       break;
		    }

		}
			    
	    }

	    if( n->lower )
	    {
	        intersectionPoint( numberIntersections, n->lower, ray, point, (dirCounter+1)%2 );
	    }

	    if( n->midSmall )
	    {
	        intersectionPoint( numberIntersections, n->midSmall, ray, point, (dirCounter+1)%2 );
	    }
	    
	    if( n->midLarge )
	    {
	        intersectionPoint( numberIntersections, n->midLarge, ray, point, (dirCounter+1)%2 );
	    }
	    
	    if( n->upper )
	    {
	        intersectionPoint( numberIntersections, n->upper, ray, point, (dirCounter+1)%2 );
	    }
	    	    
	    
	}
	
    }
    
};

#endif












