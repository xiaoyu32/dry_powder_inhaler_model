#ifndef _RAYSEARCH_
#define _RAYSEARCH_

#ifdef VECTORIZE
    #define VSIZE (4)
#else
    #define VSIZE (3)
#endif

#define RESTRICT __restrict__

#include <vector>
#include <fstream>

using namespace std;

#define NRAYS 6

#define VOTING_BIAS (0)

//intersection KD tree implementation
#include "KDtree_search.h"

class RaySearch{

   public :
   
   RaySearch();
   
   //construct from ascii file
   RaySearch( ifstream& is );
   
   ~RaySearch();	
   
   void add_face( vector<double*>&, double *);   	
   bool is_inside( const double * point ) const;	
   
   bool free_memory;
   bool print_flag;
   
   void printFaces() const;
    
   vector< vector< double* > >& get_faces(){ return faces; }
   vector< double * >& get_face_normals(){ return face_normals; }
   
   void print_faces();
   void compress(double); 
   
   void initSearchTrees();
   	
   protected :	
   
   void initObject();
   
   void getRandomPermutation( int* permutation, int N ) const;
   
   //double determinant( const double* a, const double* b , const double* c ) const;
   //void swap( double* a, double* b, int i ) const;
   //void sum( const double* a, double* b, int i ) const;
   //void gaussianSolve( double* sol, const double* a , const double* b, const double* c, const double* rhs ) const;
   //bool solveTriangle( const double*, const double*, const double*, const double*, const double* ) const;
   //int intersectionPoints( const double*, const vector< double* >&, const double* ) const;
   
   void add_face_vec( vector<double*>& face, double * normal );
   double* minus( double* a, double * b, double * c );
   
   void normalize( double* );
   
   // -- create base vectors a & b that are normal to n --
   void createBase( double* n, double* a, double* b );
      
   double norm( const double * ) const;
   double dot( double *, double * );
   
   double sign( const double& ) const; 
   double abs( const double& ) const;
   
   double angle(  double * a, double * b,
	    	  double* p1, double * p2,
		  double* p3, double * ave  );
			
   double angle2( double * a, double * b );
   double angle3( double * a, double * b, double * n, bool debug );
   
   void orientate_face( vector<double*>& face, double* n, bool& );
   
   double cross2D( double* a, double* b );
   
   bool planeCheck( vector<double*>&, vector<double*>&  );
   
   void get_normal_vector( double* p1, double* p2, double* nv );
   void cross( double * a, double * b, double * c );
   
   void clean_faces();
   void compress_face(vector<double*>& face, double* n, double tol);
   void merge( int i_, int j_ );
   void shallow_cpy( double* v1, double* v2 );
   int geti( int i, int n );
   
   void get_indexes( vector<int>& ii, vector<int>& jj, 
   		     vector<double*>& face1, vector<double*>& face2,
	             int count );
   
   int get_count(  vector<double*>& face1, vector<double*>& face2 );
   int get_count2(  vector<double*>& face1, vector<double*>& face2 );
   bool compare( double tol, vector<double*>& face1, vector<double*>& face2, double* n1, double* n2 );
   bool same_point( double* a, double* b );
   bool on_face( double* a, vector<double*>& face );
   
   vector< vector<double*> > faces;
   vector<double*> face_normals;
   
   // -- rays --
   double rays[NRAYS][VSIZE];
   
   // -- normal surfaces matching the rays --
   double rayBase[NRAYS][2][VSIZE];
   
   
   KdTreeSearch* searchTrees[NRAYS];
   
   	
};


#endif
