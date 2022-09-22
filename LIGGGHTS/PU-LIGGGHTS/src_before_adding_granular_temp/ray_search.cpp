#include "ray_search.h"
#include <vector>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>

#define small  1.0e-10
#define SMALL_ 6.07
#define PI 3.14159265359

/*
#define VECTORIZE
*/

using namespace std;

void RaySearch::initObject()
{
    free_memory = 0;
    print_flag = false;

    for( int i = 0; i < NRAYS; ++i )
        for( int j = 0; j < VSIZE; ++j )
            this->rays[i][j] = 0;
    
    for( int i = 0; i < NRAYS; ++i )
        for( int j = 0; j < 2; ++j )
	   for( int k = 0; k < VSIZE; ++k )
               rayBase[i][j][k] = 0;
    
    //this->rays[0][0] = 1;
    //this->rays[1][1] = 1;
    //this->rays[2][2] = 1;
    
    // -- ray allocations --
    this->rays[0][0] = 1;
    this->rays[0][1] = 0.04;
    this->rays[0][2] = -0.05;

    this->rays[1][0] = -0.02;
    this->rays[1][1] = 1;
    this->rays[1][2] = 0.1;
    
    this->rays[2][0] = -0.02;
    this->rays[2][1] = 0.01;
    this->rays[2][2] = 1;
    
    this->rays[3][0] = -1.;
    this->rays[3][1] = -0.08;
    this->rays[3][2] = -0.06;

    this->rays[4][0] = -0.1;
    this->rays[4][1] = -1;
    this->rays[4][2] = -0.14;

    this->rays[5][0] = -0.01;
    this->rays[5][1] = 0.02;
    this->rays[5][2] = -1;

    for( int i = 0; i < NRAYS; ++i )
       normalize( this->rays[i] );
    
    //std::cout<<"Ray00 = "<<rays[0][0]<<" "<<rays[0][1]<<" "<<rays[0][2]<<std::endl;
    //std::cout<<"Ray11 = "<<rays[1][0]<<" "<<rays[1][1]<<" "<<rays[1][2]<<std::endl;
    //std::cout<<"Ray22 = "<<rays[2][0]<<" "<<rays[2][1]<<" "<<rays[2][2]<<std::endl;  
    
    // - init search trees for each ray -
    //searchTrees[0] = new KdTreeSearch( rays[1], rays[2] ); // used for ray1
    //searchTrees[1] = new KdTreeSearch( rays[0], rays[2] ); // used for ray2
    //searchTrees[2] = new KdTreeSearch( rays[0], rays[1] ); // used for ray3
    
    
    for( int i = 0; i < NRAYS; ++i )
    {
        createBase( this->rays[i], this->rayBase[i][0], this->rayBase[i][1] );
    }
    
    
    for( int i = 0; i < NRAYS; ++i )
    {
        searchTrees[i] = new KdTreeSearch( this->rayBase[i][0], this->rayBase[i][1] );
    }
    
    
    #ifdef VECTORIZE
    cout<<"Using vectorized version."<<endl; 
    #endif
}

void RaySearch::createBase( double* n, double* a, double* b )
{
    
    double ee[3] = {0,0,0};
    
    // -- create random vector (must not be parallel to n) --
    do{
       for( int i = 0; i < 3; ++i )
       {
           ee[i] = double( rand() % 10000000 )/10000000;
       }
    }while( abs( angle2(n,ee) ) < 1e-6 );
    
    normalize( ee );
    
    cross( n, ee, a );
    cross( a, n, b );
    
    normalize( a );
    normalize( b );
    
}


void RaySearch::normalize( double* a )
{
   
   const double rr = norm(a);
   for( int i = 0; i < VSIZE; ++i )
       a[i] /= rr;
   
}

void RaySearch::getRandomPermutation( int* permutation, int N ) const
{
    
    for( int i = 0; i < N; ++i )
    {
        permutation[i] = i;
    }
    
    for( int i = 1; i < N; ++i )
    {

	int j = rand() % i;

	int ii = permutation[i];
	permutation[i] = permutation[j];
	permutation[j] = ii;

    }
    
}

void RaySearch::initSearchTrees()
{
    
    // - use random permutation to ensure that the faces are not inserted in accending or decending order -
    int* randomPermutation = new int[faces.size()];
    getRandomPermutation( randomPermutation, faces.size() );
       
    for( int i = 0; i < NRAYS; ++i )
    {
       for( int j = 0; j < faces.size(); ++j )
       {
	   searchTrees[i]->addFace( faces[ randomPermutation[j] ] );
       }
    }
    
    for( int i = 0; i < NRAYS; ++i )
       searchTrees[i]->initTree();
    
    delete[] randomPermutation;
        
}



RaySearch::RaySearch()
{
    initObject();   
}

RaySearch::RaySearch( ifstream& file )
{    

    initObject();

    string line;
    
    while( getline( file, line ) )
    {
        
	if( line.size() >= 7 && strncmp( "Normal:", line.c_str(), 7 ) == 0 )
	{
	    
	    double* normal = new double[3];
	    
	    stringstream ss;
	    ss<<line;   
	    
	    ss>>normal[0]>>normal[1]>>normal[2];
	    
	    getline( file, line );
	   
	    stringstream sv;
	    sv<<line;
	    
	    std::vector<double*> face;
	    
	    while( !sv.eof() )
	    {
	        
		double* v = new double[3];
		
		
		string temp;
		  
		for( int i = 0; i < 3; ++i )  	    
	        {  
		   sv>>temp;
		   stringstream(temp)>>v[i];
		}
		
		//std::cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<std::endl;
		
		face.push_back(v);
	    }
	    //std::cout<<std::endl;
	    add_face( face, normal );
	    
	}
	
    }
    
    file.close();
    
    //std::cout<<"Ray0 = "<<this->rays[0][0]<<" "<<this->rays[0][1]<<" "<<this->rays[0][2]<<std::endl;
    //std::cout<<"Ray1 = "<<this->rays[1][0]<<" "<<this->rays[1][1]<<" "<<this->rays[1][2]<<std::endl;
    //std::cout<<"Ray2 = "<<this->rays[2][0]<<" "<<this->rays[2][1]<<" "<<this->rays[2][2]<<std::endl;    
}

RaySearch::~RaySearch( )
{

   if( !free_memory ) return;

   //free inserted vectors
   for( int i = 0; i < faces.size(); ++i )
   {
       delete[] face_normals.at(i);
       
       for( int j = 0; j < faces.at(i).size(); ++j )
       {
           delete[] faces.at(i).at(j);
       }
   }
   
   for( int i = 0; i < NRAYS; ++i )
   {
       delete searchTrees[i];
   }
   
}

//Add face to the polyhedron
void RaySearch::add_face( vector<double*>& face, double * normal )
{
    
    if( face.size() < 3 )
    {
        cout<<"Error: Too few vertices!"<<endl;
        return;
    }
    	
    #ifndef VECTORIZE	
    
    faces.push_back( face );
    face_normals.push_back( normal );
    
    #else
    
    vector<double*> padded_face;
    double * v;
    double * padded_normal = new double[VSIZE];
    shallow_cpy( normal, padded_normal );
    
    padded_face.reserve( face.size() );
    
    for( int i = 0; i < face.size(); ++i )
    {
        
	v = new double[VSIZE];
	shallow_cpy( face.at(i), v );
	padded_face.push_back( v );
	
	if( free_memory ){
	   delete[] face.at(i);
	   face.at(i) = NULL;
	}
	
    }
    
    add_face_vec( padded_face, padded_normal );
        
    #endif
    
}

void RaySearch::add_face_vec( vector<double*>& face, double * normal )
{
    faces.push_back( face );
    face_normals.push_back( normal );  
}



double RaySearch::norm( const double * RESTRICT a ) const
{
    double res = 0.0;
    
    #ifndef VECTORIZE
    for( int i = 0; i < VSIZE; ++i )
    {
    	res += a[i]*a[i];
    }
    #else
    
    __builtin_assume_aligned(a, VSIZE);
    
    for( int i = 0; i < VSIZE; ++i )
    {
        res += a[i]*a[i];
    }
    
    #endif
    
    	
	
    return sqrt( res );	
}

double * RaySearch::minus( double* RESTRICT a, double * RESTRICT b, double * RESTRICT c )
{
    
    #ifndef VECTORIZE
    for( int i = 0; i < VSIZE; ++i )
       c[i] = a[i] - b[i];
    
    #else
    
    __builtin_assume_aligned(a, VSIZE);
    __builtin_assume_aligned(b, VSIZE);
    __builtin_assume_aligned(c, VSIZE);
    
    for( int i = 0; i < VSIZE; ++i )
       c[i] = a[i] - b[i];
    
    #endif
    
    return c;
    
}

void RaySearch::get_normal_vector( double* RESTRICT p1, double* RESTRICT p2, double* RESTRICT nv)
{
    
    cross( p1, p2, nv );
    double r = norm(nv);
    
    if( r <= 0 ){
       r = 1.0;
    }
    
    #ifdef VECTORIZE
    __builtin_assume_aligned(nv, VSIZE);
    #endif 
    
    for( int i = 0; i < VSIZE; ++i )
        nv[i] /= r;
    
}




double RaySearch::dot( double * RESTRICT a, double * RESTRICT b )
{
     
    double res = 0.0; 
    
    #ifndef VECTORIZE
    for( int i = 0; i < VSIZE; ++i )
    {
    	res += a[i]*b[i];
    }
    #else
    
    __builtin_assume_aligned(a, VSIZE);
    __builtin_assume_aligned(b, VSIZE);
    
    for( int i = 0; i < VSIZE; ++i )
    {
    	res += a[i]*b[i];
    }
    #endif	
	
	
    return res;
}

/*
double RaySearch::determinant( const double* a, const double* b , const double* c ) const
{
    
    return   a[0] * ( b[1] * c[2] - c[1] * b[2] ) 
    	   - a[1] * ( b[0] * c[2] - c[0] * b[2] )
    	   + a[2] * ( b[0] * c[1] - c[0] * b[1] );
}


void RaySearch::swap( double* a, double* b, int i ) const
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

void RaySearch::sum( const double* a, double* b, int i ) const
{
    
    if( abs( a[i] ) < small ){
    	
	//this should never happen
	std::cout<<"Gaussian Elimination Failed !!!"<<std::endl;
	
        return;
    }
    
    double cof = -b[i]/a[i];
    
    for( int k = i; k < 4; ++k )
    {
        
	b[k] += cof* a[k];
	
    }
    
    
} */

/*
void RaySearch::gaussianSolve( double* RESTRICT sol, const double* RESTRICT a , const double* RESTRICT b, const double* RESTRICT c, const double* RESTRICT rhs ) const
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
    
}*/

/*
bool RaySearch::solveTriangle( const double* RESTRICT pt, const double* RESTRICT a, const double* RESTRICT b, const double* RESTRICT p, const double* RESTRICT ray ) const
{
    
    //compute determinate of matrix [a b -ray]   
    double det = determinant( a, b, ray );
    
    //ray is parallel to the triangle plane
    if( abs( det ) < small )
    {
        
	//DO SOMETHING HERE (voting should take care of this even when point lies in the plane)
	
	return false;
    }
    
    
    //there will be intersection with the triangle plane and a line going through the point
    double sol[3];
    
    double rhs[3];
    
    for( int i = 0; i < 3; ++i )
    {
        rhs[i] = p[i] - pt[i];
    }
    
    gaussianSolve( sol, a, b, ray, rhs );
    
    //std::cout<<"Ray = "<<ray[0]<<" "<<ray[1]<<" "<<ray[2]<<std::endl;
    //std::cout<<"Sol = "<<sol[0]<<" "<<sol[1]<<" "<<sol[2]<<std::endl;
    //std::cout<<std::endl;
    
    //if the intersection point lies within the triangle
    if( sol[2] < small && sol[0] > -small && sol[1] > -small  && (sol[0] + sol[1] <= 1.0 ) )
    {
        return true;
    }else{
        return false;
    }
    
} */

/*
int RaySearch::intersectionPoints( const double* RESTRICT point, const vector< double* >& face, const double* RESTRICT ray ) const 
{
    
    //std::cout<<"Ray = "<<ray[0]<<" "<<ray[1]<<" "<<ray[2]<<std::endl;
    
    if( face.size() < 3 ) return 0;
    
    const double* p1 = face[0];
    double p2[VSIZE];
    double p3[VSIZE];
        
    for( int i = 2; i < face.size(); ++i )
    {
        
	for( int k = 0; k < VSIZE; ++k )
	{
	    p2[k] = face[i-1][k] - p1[k];
	    p3[k] = face[i][k] - p1[k];
	}
	
	if( solveTriangle( p1, p2, p3, point, ray ) )
	{
	   return 1;
	}
	
    }
    
    
    
    return 0;
    
}*/

bool RaySearch::is_inside( const double * point ) const
{
    
	
     int nintersections[NRAYS] = {0,0,0};
     
     for( int i = 0; i < NRAYS; ++i )
     {
         /*
         for( int j = 0; j < faces.size(); ++j )
	 {
	     nintersections[i] += intersectionPoints( point, faces.at(j), rays[i] );
	 }
	 */
	 nintersections[i] = searchTrees[i]->intersectionPoints( rays[i], point ); 
	 //std::cout<<std::endl;
     }
     
     
     //std::cout<<"Number of faces: "<<faces.size()<<std::endl;
     
     //votes (>0 means inside and <0 means outside), each ray votes once
     int inside = 0;
          	
     for( int i = 0; i < NRAYS; ++i )	
     {
         
	 //std::cout<<"Intersections: "<<i<<" "<<nintersections[i]<<std::endl;
	 
	 if( nintersections[i]%2 == 0 )
	 {
	     --inside;
	 }else{
	     ++inside;
	 } 
	 
     }	
     
     //std::cout<<"Outcome: "<<inside<<std::endl;
     
     // -- bias voting towards positive votes --	
     return inside>=VOTING_BIAS;	

}

double RaySearch::cross2D( double* a, double* b )
{
    return a[0]*b[1] - a[1]*b[0];
}


void RaySearch::cross( double * RESTRICT a, double * RESTRICT b, double * RESTRICT c )
{
        
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    
    #ifdef VECTORIZE
    c[3] = 0.0;
    #endif
    
}


//hejssan
double RaySearch::angle( double * RESTRICT a , double * RESTRICT b, 
				double* RESTRICT p1, double * RESTRICT p2, 
				double* RESTRICT p3, double * RESTRICT ave  )
{

    double res = dot(a,b);
    
    const double r1 = norm(a);
    const double r2 = norm(b);
    
    if( r1 <= 0 || r2 <= 0 ){
	return 0;
    }
    
    res /= r1; 
    res /= r2;
    
    if( res > 1.0 ) res = 1.0;
    if( res < -1.0 ) res = -1.0;
        
    res = acos( res );
    
    //actual angle between the planes
    res = PI - res;
   
    //check if the corner is concave
    double aa[VSIZE];
    double bb[VSIZE];
    double cc[VSIZE];
    
    #ifdef VECTORIZE
    __builtin_assume_aligned(p1, VSIZE);
    __builtin_assume_aligned(p2, VSIZE);
    __builtin_assume_aligned(p3, VSIZE);
    #endif
    
    for( int i = 0; i < VSIZE; ++i )
    {
        aa[i] = p2[i] - p1[i];
	bb[i] = p3[i] - p2[i];
    }
    
    cross( aa, bb, cc );
    
    double res2 = angle2(aa,bb);
    
    //the non projected corner is almost straight
    if( abs(res2) < 0.01 )
    {
        return PI;
    }
    
    //r1 = norm(cc);
    //for( int i = 0; i < 3; ++i )
    //    cc[i] /= r1;
    
    double ss = dot( ave, cc );
    
    //concave corner
    if( ss < 0  ){
       res = 2*PI-res; 	
    }
    	
    //if( res < 0 ) res = 2*PI-res;
    //while( res > 2*PI ) res = res-2*PI;
        	
    return res;
}




double RaySearch::angle2( double * a, double * b )
{

    
    double res = dot(a,b);
    
    double r1 = norm(a);
    double r2 = norm(b);
    
    if( r1 <= 0 || r2 <= 0 ){
	return 0;
    }
    
    res /= r1; 
    res /= r2;
    
    if( res > 1.0 ) res = 1.0;
    
    if( res < -1.0 ) res = -1.0;
    
    res = acos( res );
    	
    return abs(res);
}

double RaySearch::angle3( double * a, double * b, double* n, bool debug = false )
{

    
    double res = dot(a,b);
    
    double r1 = norm(a);
    double r2 = norm(b);
    
    if( r1 <= 0 || r2 <= 0 ){
	return 0;
    }
    
    res /= r1; 
    res /= r2;
    
    if( res > 1.0 ) res = 1.0;
    
    if( res < -1.0 ) res = -1.0;
    
    res = acos( res );
    
    //interior angle
    res = PI-res;
    
        
    double cc[VSIZE];
        
    cross( a, b, cc );	
    
    double ss = dot( cc, n );
    
    //concave corner
    if( ss < 0 ){
       //if( print_flag && debug ) cout<< "Concave corner detected! "<<ss<<endl;
       res = 2*PI-res; 	
    }	
    
    if( res < 0 ) res = 2*PI - res;
    while( res > 2*PI ) res = res - 2*PI;
    
    //if( print_flag ) cout<< "Angle "<< res <<endl;
        
    return res;
}

double RaySearch::sign( const double& a ) const 
{
    if( a >= 0 ) return 1.0;
    else	 return -1.0;
}

double RaySearch::abs( const double& a ) const
{
    return sign(a)*a;
}


void RaySearch::printFaces() const
{
    
    for( int i = 0; i < faces.size(); ++i )
    {
    
        std::cout<<"Face "<<i<<":"<<std::endl;
        for( int j = 0; j < faces[i].size(); ++j )
	{
	    std::cout   <<"( "<<faces[i][j][0]<<" "
	    		<<faces[i][j][1]<<" "
			<<faces[i][j][2]<<" ) ";
	
	}
    	std::cout<<std::endl;
    
    }
    
    
}




//FACE COMPRESSION 

void RaySearch::compress( double tol )
{
    
    bool cflag = true;
    
    vector< bool > merge_list;
    
    int n_orig_ = faces.size();
    int n_orig;
    
    int it = 0;
    int max_it = 200;
    
    while( cflag && it < max_it )
    {
       
       cflag = false;
       n_orig = faces.size();
       
       
       merge_list.clear();
       merge_list.resize( n_orig, false );
       
       for( int i = 0; i < n_orig; ++i )
       {
    	   	   
	   if( merge_list.at(i) ) continue;	   
		   
	   for( int j = i; j < n_orig; ++j )
	   {
		
		if( i==j ) continue;
		
		if( merge_list.at(j) ) continue;
		
		if( compare( tol, faces.at(i), faces.at(j), 
		             face_normals.at(i), face_normals.at(j) ) 
			   )
		{
		
		    cflag = true;
		    merge( i, j );
		    merge_list.at(i) = true;
		    merge_list.at(j) = true;
		    break;
		    
		}
		   
    	   }
	   
       }
       
       
       ++it;
       clean_faces();
       //if( print_flag )  cout<<"Current number of faces: "<<faces.size()<< " " <<cflag<<endl;
       
    }
    
    if( it == max_it -1 && print_flag ) cout<<"Warning: Iteration limit reached!"<<endl;
    
    //clean_faces();
    for( int i = 0; i < faces.size(); ++i )
        compress_face( faces.at(i), face_normals.at(i), tol );
    
    
    if( print_flag ) cout<<"Compressed number of faces from "<<n_orig_<<" to "<< faces.size()<< endl;
    
    
    //insert compressed faces into the search trees
    //initSearchTrees();
    
}

void RaySearch::clean_faces()
{
    
    int nn = faces.size();
    int flag = true;
    
    
    while( flag ){
    	
	flag = false;
	
    	for( int i = 0; i < faces.size(); ++i )
    	{
        
	    if( faces.at(i).size() == 0 )
	    {
	        flag = true;
		faces.erase( faces.begin()+i );
		delete[] face_normals.at(i);
		face_normals.erase( face_normals.begin()+i );
            }
	
    	}
	
    }
          
}

void RaySearch::compress_face(vector<double*>& face, double* n, double tol)
{
    
    bool cflag = true;
    int n_orig = face.size();
    
    double p1[VSIZE];
    double p2[VSIZE];
    
    
    if( face.size() == 0 ) return;
    
    while( cflag )
    {
    
        cflag = false;
        
	for( int i = 1; i < face.size()-1; ++i )
	{
	       
	    minus( face.at(i), face.at(i-1), p1 );
	    minus( face.at(i+1), face.at(i), p2 );
	    
	    if( angle2(p1, p2 ) < tol ){
	        cflag = true;
	        
            if( free_memory ){
	        if( face.at(i) ) delete[] face.at(i);
	    }
		
		face.erase( face.begin() + i );
		
		if( face.size() < 3 ) cout<<"Error: tolerance is too rough!"<<endl;
	    }
    
	}
	
	//first and last vertex
	minus( face.at( face.size()-1 ), face.at( face.size()-2 ), p1 );
	minus( face.at(0), face.at( face.size()-1 ), p2 );    

      	if( angle2(p1, p2 ) < tol ){
	    cflag = true;
	    
            if( free_memory ){
	        if( face.at(face.size()-1) ) delete[] face.at(face.size()-1);
	    }
	    
	    face.erase( face.begin() + face.size()-1 );

	    if( face.size() < 3 ) cout<<"Error: tolerance is too rough!"<<endl;
	}
	
	minus( face.at(0), face.at( face.size()-1 ), p1 );
	minus( face.at(1), face.at( 0 ), p2 );    

      	if( angle2(p1, p2 ) < tol ){
	    cflag = true;
	    
	    if( free_memory ){
	        if( face.at(0) ) delete[] face.at(0);
	    }
	    
	    face.erase( face.begin() );

	    if( face.size() < 3 ) cout<<"Error: tolerance is too rough!"<<endl;
	}
	
    }
    
    bool success = true;
    
    orientate_face( face,  n, success );
          
}

bool RaySearch::planeCheck( vector<double*>& face1, vector<double*>& face2 )
{
    
    if( face1.size() < 3 || face2.size() < 3 ) return false;
    
    double* p1 = face1.at(0);
    
    double p2[3] = {0,0,0};
    double p3[3] = {0,0,0};
    
    minus( p1, face1.at(1), p2 );
    minus( p1, face1.at(2), p3 );

    double normal[3] = {0,0,0};
    
    cross( p1,p2,normal );
    
    double rr = norm( normal );
    
    for( int i = 0; i < 3; ++i )
        normal[i] /= rr;
    
    double distExtra = 0;

    if( face1.size() >= 4 )
    {
        minus( p1, face1.at(3), p2 );
        distExtra = dot( normal, p2 );
	distExtra = distExtra>0?distExtra:-distExtra;
    }
    
    for( int i = 0; i < face2.size(); ++i )
    {
        minus( p1, face2.at(i), p2 );
	double d = dot( normal , p2 );
	
	if( (d>0?d:-d) >= small + distExtra ) return false;
    }

    return true;
    
}

void RaySearch::merge( int i_, int j_ )
{
       
    if( faces.at(i_).size() == 0 || faces.at(j_).size() == 0 ){
        cout<<"Error: nothing to merge!"<<endl;
	return;
    } 
    
    int i,j;
    vector<int> ii;
    vector<int> jj;
    
    int si = faces.at(i_).size();
    int sj = faces.at(j_).size();
    
    
    ii.reserve( si );
    jj.reserve( sj );
    
    int count = get_count( faces.at(i_), faces.at(j_) );
    
    if( count < 2 ) return;
    
    
    if( !planeCheck( faces.at(i_), faces.at(j_ ) ) ) return;
    
    
    vector< double* > new_face;
    double* v;
    
    double* n1 = face_normals.at(i_);
    double* n2 = face_normals.at(j_);
    
    if( count >= 2 )
    {
        
	if( faces.at(i_).size() == count ){
	
	    cout<<"Error"<<endl;
	    if( free_memory ){

		for( i = 0; i < si; ++i )
		{
		    if( faces.at(i_).at(i) ) delete[] faces.at(i_).at(i);
		    faces.at(i_).at(i) = NULL;
		}

	    }
	    
	    faces.at(i_).clear();
	    
	    return;
	    
	}
	
	if( faces.at(j_).size() == count ){
 	    
	    cout<<"Error"<<endl;
	    
	    if( free_memory ){

		for( j = 0; j < sj; ++j )
		{
		    if( faces.at(j_).at(j) ) delete[] faces.at(j_).at(j);
		    faces.at(j_).at(j) = NULL;
		}

	    }
	    
	    faces.at(j_).clear();
	    
	    return;
	}
	
	
        new_face.reserve( faces.at(i_).size() + faces.at(j_).size() + 1 - count );
        get_indexes( ii, jj, faces.at(i_), faces.at(j_), count );
	
	if( ii.size() == 0 || jj.size() == 0 )
	{
	    cout<< "Error: no vertices found!"<<endl;
	    return;
	}
	
	
	for( i = 0; i < ii.size(); ++i )
	{
	   v = new double[VSIZE];
	   shallow_cpy( faces.at(i_).at( ii.at(i) ), v );
	   new_face.push_back( v );
	}
	
	for( j = 0; j < jj.size(); ++j )
	{
	   v = new double[VSIZE];
	   shallow_cpy( faces.at(j_).at( jj.at(j) ), v );
	   new_face.push_back( v );
	}
	
	//cout<< new_face.size()<< " " << count<< " " << si << " " << sj << endl; 

	//TODO: implement face normal recomputation
	v = new double[VSIZE];
	shallow_cpy( n2, v );
	
	bool success_flag = false;
	//ensure that the new face is oriented clockwise
	orientate_face( new_face, v, success_flag );
	
	//add new face
	add_face_vec( new_face, v );
		
	if( free_memory ){
	    
	    for( i = 0; i < si; ++i )
	    {
	        if( faces.at(i_).at(i) ) delete[] faces.at(i_).at(i);
		faces.at(i_).at(i) = NULL;
	    }
	    
	    for( j = 0; j < sj; ++j )
	    {
	        if( faces.at(j_).at(j) ) delete[] faces.at(j_).at(j);
		faces.at(j_).at(j) = NULL;
	    }
	    
	}
	
        faces.at(i_).clear();
	faces.at(j_).clear();

	
    }else{
        cout<<"Error: vertex count incompatible!"<<endl;
    }
    
    
}

inline void RaySearch::shallow_cpy( double* v1, double* v2 )
{
    
    if( !v1 ) 
    {
        cout<<"Error: Null pointer!"<<endl;
	return;
    }

    v2[0] = v1[0];
    v2[1] = v1[1];
    v2[2] = v1[2];
    
    #ifdef VECTORIZE
    v2[3] = 0;
    #endif
    
}

inline int RaySearch::geti( int i, int n )
{
  
    if( i < 0 ) return n + i;
    if( i >= n ) return i - n;
    
    return i;
}

void RaySearch::get_indexes( vector<int>& ii, vector<int>& jj, 
				vector<double*>& face1, vector<double*>& face2, 
				int count )
{
    
    double buff[VSIZE];
    
    double icut[2][VSIZE];	
    double jcut[2][VSIZE];	
    	    
    int ic[2];
    int jc[2];
    
    int i,j,k;
    bool flag[2];
    
    for( k = 0; k < 2; ++k )
    {
        ic[k] = -1;
	jc[k] = -1;
	flag[k] = false;
    }
        
    int si,sj;
    int ind0, ind1;
    
    si = face1.size();
    sj = face2.size();
    
    int count2 = get_count( face1, face2 );

    for( i = 0; i < si; ++i )
    {
	
	if( on_face( face1.at(i), face2 ) )
	{
	    
	    if( !on_face( face1.at( geti(i-1,si) ), face2 ) && !flag[0] ){
	        
		shallow_cpy( face1.at(i), icut[0] );
		ic[0] = i;
		flag[0] = true;
		
	    }else if( !on_face( face1.at( geti(i+1,si) ), face2 ) && !flag[1] ){
	        shallow_cpy( face1.at(i), icut[1] );
		ic[1] = i;
		flag[1] = true;
	    }
	    
	    
	}
	
	if( flag[0] && flag[1] ) break;
	
    }
    
    count2 = get_count( face1, face2 );

    if( !flag[0] || !flag[1] ){ cout<<"Error!"<<endl;  return; }

    flag[0] = false;
    flag[1] = false;
    
    
    count2 = 0;
    
    for( i = 0; i < sj; ++i )
    {
	
	if( on_face( face2.at(i), face1 ) )
	{
	    ++count2;
	    if( !on_face( face2.at( geti(i-1,sj) ), face1 ) && !flag[0] ){
	        
		shallow_cpy( face2.at(i), jcut[0] );
		jc[0] = i;
		flag[0] = true;
		
	    }else if( !on_face( face2.at( geti(i+1,sj) ), face1 ) && !flag[1] ){
	        shallow_cpy( face2.at(i), jcut[1] );
		jc[1] = i;
		flag[1] = true;
	    }
	    
	    
	}
	
	if( flag[0] && flag[1] ) break;
	
    }
    
    count2 = get_count( face1, face2 );

    if( !flag[0] || !flag[1] ){ 
    
        cout<<"Error2!"<<flag[0]<<flag[1]<< " " << count<<" "<< count2 <<endl;
	
	count = get_count( face1, face2 );
	
	cout<< "COUNT : "<< count << endl;
	
        return;
    }
    
    i = ic[1];
    while( i != ic[0] )
    {
        ii.push_back(i);
	i = geti( i+1, si );
    }
    
    int step = 0;
    
    if( same_point( face1.at(i) , jcut[0] ) )
    {
        ind0 = jc[0];
	ind1 = jc[1];
	step = -1;
    }else if( same_point( face1.at(i) , jcut[1] ) ){
        ind0 = jc[1];
	ind1 = jc[0];
	step = 1;
    }else{
        cout<<"Error3!!!!!!!!!!<<endl";
	return;
    }
    
    j = ind0;
    while( j != ind1 )
    {
         jj.push_back(j);    
     	 j = geti( j + step, sj );
    }
     
    
    
}

bool RaySearch::on_face( double* a, vector<double*>& face )
{
    
    for( int i = 0; i < face.size(); ++i )
    {
        
	if( same_point( a, face.at(i) ) )return true;
	
    }
    
    return false;
    
}

bool RaySearch::same_point( double* a, double* b )
{
    double buff[VSIZE];
    
    if( norm( minus( a, b, buff ) ) < small )
    {
        return true;   
    }
    
    return false;
}

int RaySearch::get_count2( vector<double*>& face1, vector<double*>& face2 )
{
    int i, j, count;
    double buff[VSIZE];
    count = 0;
    
    for( i = 0; i < face1.size(); ++i )
    {
        
	for( j = 0; j < face2.size(); ++j )
	{

	    if( norm( minus( face1.at(i), face2.at(j), buff ) ) < small )
	    {
	        ++count;
		break;
	    }

	}
	
    }
    
    return count;
}

int RaySearch::get_count(  vector<double*>& face1, vector<double*>& face2 )
{
    
    int count = 0;
    double buff[VSIZE];
    bool flag = false;
    bool flag2 = false;
    
    int i,j;
    
    for( i = 0; i < face1.size(); ++i )
    {
    
        flag2 = false;
        for( j = 0; j < face2.size(); ++j )
	{

	    if( norm( minus( face1.at(i), face2.at(j), buff ) ) < small )
	    {
	        flag2 = true;
	        flag = true;
	        ++count;
		break;
	    }

	}
	
	if( flag && !flag2 ){
	    
	    if( i - count > 0 ) return count;
	    
	    for( i = face1.size()-1; i > 0; --i ){
	        flag2 = false;
		
	        for( j = 0; j < face2.size(); ++j ){
		    if( norm( minus( face1.at(i), face2.at(j), buff ) ) < small )
		    {
	        	flag2 = true;
	        	++count;
			break;
		    }
		}
		
		if( flag && !flag2 ) return count;
	    }
	    
	    return count;
	    
	}
	
    }
    
    return count;
    
}


bool RaySearch::compare( double tol, vector<double*>& face1, vector<double*>& face2, double* n1, double* n2 )
{
            
    //double angle_ = angle2( n1, n2 );   
       
    //if( angle_ > tol && abs(angle_-PI) > tol )  return false;

    if( !same_point( n1, n2 ) ) return false;

    int count = get_count( face1, face2 );
    int count2 = get_count2( face1, face2 );
       
    if( count != count2 ) return false;   
       
    if( count >= 2 ) return true;
    
    return false;
    
}


void RaySearch::print_faces()
{
    
    cout<< "Face information: "<<endl;
    
    for( int i = 0; i < faces.size(); ++i )
    {
        cout<<"Face: "<<faces.at(i).size()<<" : ";
	for( int j = 0; j < faces.at(i).size(); ++j )
	{
	    //if( abs( abs( face_normals.at(i)[2]) -1 ) < small  ){
	    cout<< faces.at(i).at(j)[0] << ", ";
	    cout<< faces.at(i).at(j)[1] << ", ";
	    cout<< faces.at(i).at(j)[2] << "  ";
	    cout<< " ,";
	    //}
	}
	
	
	cout<<"Normal : "<< face_normals.at(i)[0] << ", "<< face_normals.at(i)[1]<<", "<< face_normals.at(i)[2]<<endl;
	
	cout<<endl;
    }
    
}


void RaySearch::orientate_face( vector<double*>& face, double* n, bool& success )
{
    
    //compute interior angle sum of the non projected polygon
    double interior = 0;
    int l = face.size();
    
    double a[3];
    double b[3];
    double angle;
    vector<double*> vbuffer;
    
    for( int k = 0; k < 3; ++k )
    {
        a[k] = face.at(0)[k] - face.at(l-1)[k];
	b[k] = face.at(1)[k] - face.at(0)[k];
    }
    
    angle = angle3( a,b,n, success ); 
    interior += angle;
    
    for( int i = 1; i < l-1; ++i )
    {
	
	for( int k = 0; k < 3; ++k )
	{
	    a[k] = face.at(i)[k] - face.at(i-1)[k];
	    b[k] = face.at(i+1)[k] - face.at(i)[k];
	}
	
	angle = angle3( a,b,n, success ); 
        interior += angle; 
	   
    }
    
    for( int k = 0; k < 3; ++k )
    {
        a[k] = face.at(l-1)[k] - face.at(l-2)[k];
	b[k] = face.at(0)[k] - face.at(l-1)[k];
    }
    
    angle = angle3( a,b,n, success ); 
    interior += angle;
    
    /*if( print_flag && success){
	cout<<" "<<endl;
	cout<<interior<<endl;
	cout<<(l-2)*PI<<endl;
	cout<<(l+2)*PI<<endl;
	cout<<"Angles!"<<endl;
    }*/ 
    
    //reorient the face
    if( abs( interior - (l-2)*PI ) > SMALL_ && abs( interior - (l+2)*PI ) < SMALL_  )
    {

	vbuffer.reserve( l );
	
	for( int i = 0; i < l; ++i )
	{
	    vbuffer.push_back( face.at(l-1-i) );
	}
	
	face = vbuffer; 
	
	if( print_flag && success){
	    cout<<" "<<endl;
	    cout<<interior<<endl;
	    cout<<(l-2)*PI<<endl;
	    cout<<(l+2)*PI<<endl;
	    cout<<"Orientation changed!"<<endl;
	} 
	
	success = true;
	
    }else if( abs( interior - (l-2)*PI ) < SMALL_ ){
        
	success = true;
	
    }else{
        
	if( print_flag && success){
	    cout<<" "<<endl;
	    cout<<interior<<endl;
	    cout<<(l-2)*PI<<endl;
	    cout<<(l+2)*PI<<endl;
	    cout<<"Failed Merging!"<<endl;
	} 
	
        success = false;
    
    }  
    /*
    if( !success ) return;
    
    vbuffer.reserve( l );

    for( int i = 0; i < l; ++i )
    {
	vbuffer.push_back( face.at(l-1-i) );
    }
    
    for( int k = 0; k < 3; ++k )
    {
        a[k] = vbuffer.at(0)[k] - vbuffer.at(l-1)[k];
	b[k] = vbuffer.at(1)[k] - vbuffer.at(0)[k];
    }
    
    angle = angle3( a,b,n, success ); 
    interior = angle;
    
    for( int i = 1; i < l-1; ++i )
    {
	
	for( int k = 0; k < 3; ++k )
	{
	    a[k] = vbuffer.at(i)[k] - vbuffer.at(i-1)[k];
	    b[k] = vbuffer.at(i+1)[k] - vbuffer.at(i)[k];
	}
	
	angle = angle3( a,b,n, success ); 
        interior += angle; 
	   
    }
    
    for( int k = 0; k < 3; ++k )
    {
        a[k] = vbuffer.at(l-1)[k] - vbuffer.at(l-2)[k];
	b[k] = vbuffer.at(0)[k] - vbuffer.at(l-1)[k];
    }
    
    angle = angle3( a,b,n, success ); 
    interior += angle;
    
    if( print_flag && success){
	cout<<" "<<endl;
	cout<<interior<<endl;
	cout<<(l-2)*PI<<endl;
	cout<<(l+2)*PI<<endl;
	cout<<"Angles Opposite!"<<endl;
    } 
    */
    
}























