#ifndef WALL_DATA_EXCHANGEMODEL
#define WALL_DATA_EXCHANGEMODEL



#include "math.h"
#include "fvCFD.H"
#include "cfdemCloud.H"
#include <string.h>
#include <utility>
#include "boundaryElementFinder.H"




namespace Foam
{


// -- dem wall data container --
class WallValueContainer
{
   
   protected :
   
     // -- number of meshes --
     int nmesh;

     // -- number of elements in each mesh -- 
     int* nelements;

     // -- length of vector --
     int len;

     // -- number of vectors --
     int nvec;

     double**** data;
     
     char* id;
     
     int index;
     int indexMesh;
     
     void deleteData()
     {

	if( !data ) return; 

	for( int i = 0; i < nmesh; ++i )
	{
	   deleteMesh( i );
	}

	delete[] data;

	nmesh = 0;
	delete[] nelements;
	nelements = NULL;
	data = NULL;
	nvec = 0;
	len = 0;

     } 

     void deleteMesh( int imesh )
     {
		
	if( !data ) return;

	double*** meshData = data[imesh]; 

	if( !meshData ) return;

	if( nelements[imesh] == 0 ) return;

	for( int j = 0; j < nelements[imesh]; ++j )
	{

	   for( int k = 0; k < nvec; ++k )
	   {
	      delete[] data[imesh][j][k];
	   }

	   delete[] data[imesh][j];

	}      
	
	delete meshData;
	data[imesh] = NULL;
	
	nelements[imesh] = 0;
		
     }
    
   
   public :
   
     WallValueContainer( char* id_ ) :
     data( NULL ),
     nmesh( 0 ),
     nelements( NULL ),
     len( 0 ),
     nvec( 0 ),
     index( 0 ),
     indexMesh( 0 ),
     id( NULL )
     {
         if( id_ )
	 {
	     id = new char[ strlen( id_ )+1 ];
	     strcpy( id, id_ );
	 }
     }

     ~WallValueContainer()
     { 
	if( id ) delete[] id;	
	if( data ) deleteData();	 
     }

     void allocateData( int nmesh_ )
     {
	
	// -- correct number of meshes already allocated --
	if( nmesh_ == nmesh ) return;
	
	if( data ) deleteData();
	
	nmesh = nmesh_;
	
	data = new double***[nmesh];
	nelements = new int[nmesh];
	
	
	for( int i = 0; i < nmesh; ++i )
           nelements[i] = 0;

     }
     
     void begin()
     {
        index = 0;
	indexMesh = 0;
     }

     double** nextElement(  int& imesh, int& ielement )
     {
     
        // -- there are no meshes --
        if( nmesh == 0 ) return NULL;
     
        while( index >= nelements[indexMesh] )
	{
	   index = 0;
	   ++indexMesh;
	   
	   if( indexMesh >= nmesh ) return NULL;
	   
	}
        
	if( indexMesh >= nmesh ) return NULL;
        
	imesh = indexMesh;
	ielement = index;
	
        return data[indexMesh][index++];
     }

     double** nextElement()
     {
        // -- there are no meshes --
        if( nmesh == 0 ) return NULL;     
     
        while( index >= nelements[indexMesh] )
	{
	   index = 0;
	   ++indexMesh;
	   
	   if( indexMesh >= nmesh ) return NULL;
	   
	}
        
	if( indexMesh >= nmesh ) return NULL;
	
        return data[indexMesh][index++];
     }

     inline double* nextVector(  int& imesh, int& ielement )
     {
	double** buf = nextElement(imesh, ielement );
        return buf ? buf[0] : NULL;
     }
     
     inline double* nextVector()
     {
	double** buf = nextElement();
	return buf ? buf[0] : NULL;	
     }
     
     inline double nextScalar()
     {   
	double** buf = nextElement();
	return buf ? buf[0][0] : 0;  
     }
     
     int size() const
     {
         
	 int nn = 0;
	 
	 for( int i = 0; i < nmesh; ++i )
	 {
	     nn += nelements[i];
	 }
	 
	 return nn;
	 
     }

     void allocateMesh( int imesh, int nelements_, int len_, int nvec_ )
     {

	if( !data ) return;
		
	if( nelements[imesh] !=  nelements_ || len != len_ || nvec != nvec_ )
	{
            deleteMesh( imesh );
	    
	    nelements[imesh] = nelements_;
	    len = len_;
	    nvec = nvec_;

	    data[imesh] = new double**[nelements_];

	    double*** meshData = data[imesh];

	    for( int i = 0; i < nelements[imesh]; ++i )
	    {

		meshData[i] = new double*[nvec];

		for( int j = 0; j < nvec; ++j )
	           meshData[i][j] = new double[len];

	    }

	}

     }
     
     inline int getLen() const
     {
         return len;
     }
     
     inline int getNvec() const
     {
         return nvec;
     }
     
     inline char* getID()
     {
         return id;
     }
     
     inline double**** ref()
     {
	 return data;
     }
   
};



// -- wall data exchange model --
class wallDataExchangeModel
{

   public :
   
     wallDataExchangeModel( void* lmp_, cfdemCloud& sm );

     ~wallDataExchangeModel();

    // Selector

     static Foam::autoPtr<wallDataExchangeModel> New
     (
         void* lmp_, 
	 Foam::cfdemCloud& sm 
     );
     
     Foam::scalar computeArea( const double** face, int nvec ) const;
     
     void update();
     
     // -- manually initialize container --
     void initData( WallValueContainer& container );
     
     void getData( WallValueContainer& container );
     void giveData( WallValueContainer& container );
     
     // -- methods for transferring data from DEM container to OpenFoam object --
     void transferDataToField( WallValueContainer& container, surfaceScalarField& field );
     void transferDataToField( WallValueContainer& container, volScalarField& field );
     
     // -- methods for transferring data from OpenFoam object to DEM container --
     void transferDataFromField( WallValueContainer& container, surfaceScalarField& field );
     void transferDataFromField( WallValueContainer& container, volScalarField& field );
     
   protected :
     
     boundaryElementFinder* finderCFD;
     
     cfdemCloud& particleCloud_;
     void* lmp;
     
     WallValueContainer nodeCenters;
     WallValueContainer faces;
     
};

}; // end of namespace




#endif














