#ifndef _BOUNDARY_ELEMENT_FINDER_H_
#define _BOUNDARY_ELEMENT_FINDER_H_

#include <vector>
#include "fvCFD.H"
#include <cstdio>

#define INF 	(1e9)

class boundaryElementFinder
{

   protected :
    
    class Face
    {
        
	public :
	
	 Face( Foam::label patchI_, Foam::label faceI_ );
	 
	 Face();
	 
	 ~Face();
	 
	 void setIDs( Foam::label patchI_, Foam::label faceI_ );
	 
	 // -- check if point is on the face --
	 bool isInside( const Foam::vector& point, bool debug = false ) const;
	
	 void addReference( int imesh, int ielement );
	
	 void clearReferences();
	
	 void getBounds( double* min, double* max ) const;
	
	 // -- inline functions --
	 	
	 inline void addPoint( Foam::vector point )
	 {
	     points.push_back( point );   
	 }

	 inline Foam::label getPatchID() const
	 {
	     return patchI;
	 }

	 inline Foam::label getFaceID() const
	 {
	     return faceI;
	 }
	 
	 inline std::vector<int>* getMeshList()
	 {
	     return &(this->meshList);
	 }

	 inline std::vector<int>* getElementList()
	 {
	     return &(this->elementList);
	 }
	 
	 inline int size() const
	 {
	     return points.size();
	 }
	 
	 inline int nreferences() const
	 {
	     return elementList.size();
	 }  
	 
	protected :
	 
	 // -- compute counter clockwise angle of vectors a and b --
	 Foam::scalar vectorAngle( const Foam::vector& a, const Foam::vector& b ) const;
	 
	 inline Foam::scalar abs( const Foam::scalar& a ) const
	 {
	     return a > 0 ? a : -a;
	 }
	 
	 // mesh and element indices of the DEM side 
	 std::vector<int> meshList;
	 std::vector<int> elementList;
	 
	 std::vector<Foam::vector> points;
	 Foam::label patchI; 
	 Foam::label faceI;
	 
	 
	 // -- maximum allowed distance from the plane spanned by the faces (should be ASCII %f accuracy) --
	 static const double distanceThreshold = 1e-5;
	 
	 // -- maximum allowed difference from the full angle (in radians) --
	 static const double angleThreshold = 1e-2;
	 
	 // -- full angle in radians --
	 static const double fullAngle = 2 * 3.1415926;
	 
    };
    
    class Container
    {
        
	public :

	 Container( double* lower_, double* upper_ );
	 
	 ~Container();
	 
	 bool findFace( const Foam::vector& point, Foam::label& faceId, Foam::label& patchId );
	 
	 bool addReference( const Foam::vector& point, int imesh, int ielement );
	 
	 // -- returns list of all the DEM side mesh elements that belong to the fluid mesh element containing the "point" --
	 bool getDEMreferences( const Foam::vector& point, std::vector<int>*& meshList, std::vector<int>*& elementList );
	 
	 
	 // -- add face in the containers -- return fase if face does not fit in the container
	 bool addFace( boundaryElementFinder::Face& ff );
	 
	 void clearAll();
	 
	 void clearReferences();
	 
	 int countEmptyFaces() const;
	 
	 void getFace( const Foam::vector& point, Foam::label patchI, Foam::label faceI );
	 
	protected :
	 
	 
	 void addNewFace( boundaryElementFinder::Face& ff );
	 
	 void set( double* lower_, double* upper_ );

	 // -- check if the face is fully enclosed in the container --
	 bool isContained( const boundaryElementFinder::Face& face ) const;
	 
	 // -- return index of the octoTree member, which contains face -- returns -1 if not contained in any of the sub-containers
	 int containedIndex( const boundaryElementFinder::Face& face ) const;
	 
	 // -- return index of the sub-container, -1 if point is not in the subcontainer
	 int pointIndex( const Foam::vector& point ) const;
	 
	 // -- container bounds --
	 double lower[3];
	 double upper[3];
	 double mid[3];
	
	 // -- list of faces in the container that do not fit in any of the sub containers
	 std::vector< boundaryElementFinder::Face* > faces;
	 
	 // -- pointers to sub-containers --
	 Container* octoTree[8];
	
    };
    
    int nelements;
    
    boundaryElementFinder::Face** faces; 
    Foam::scalar** demFaceAreas;
    
    int npatch_;
    int* nfaces_;
    
    
    Container treeSearch;
    
   public :

    boundaryElementFinder( double* lower_, double* upper_ );
    
    void initPatches( int npatch );
    void initPatch( int ipatch, int nfaces );
    
    ~boundaryElementFinder();
    
    inline int countEmptyFaces() const
    {
        return treeSearch.countEmptyFaces();
    }
    
    // OpenFoam faces
    bool addFace( const Foam::face& face_, const Foam::pointField& face_points, Foam::label patchI, Foam::label faceI );
    
    inline int npatches() const 
    {
       return npatch_;
    }

    inline int nfaces(int i) const 
    { 
       if( i < 0 || i >= npatch_ ) return 0;
       return nfaces_[i];
    }
    
    inline double& faceArea( Foam::label patchI, Foam::label faceI )
    {
       return demFaceAreas[patchI][faceI];
    }
    
    inline bool findFace( const Foam::vector& point, Foam::label& faceId, Foam::label& patchId )
    {
       return treeSearch.findFace( point, faceId, patchId );
    }
    
    inline bool addReference( const Foam::vector& point, int imesh, int ielement )
    {
        return treeSearch.addReference( point, imesh, ielement );
    }
    
    inline void getFace( const Foam::vector& point, Foam::label patchI, Foam::label faceI )
    {
        treeSearch.getFace( point, patchI, faceI );
    }
    
    //inline bool getDEMreferences( const Foam::vector& point, std::vector<int>*& meshList, std::vector<int>*& elementList )
    //{
    //    return treeSearch.getDEMreferences( point, meshList, elementList );
    //}
    
    bool getDEMreferences( Foam::label patchI, Foam::label faceI, std::vector<int>*& meshList, std::vector<int>*& elementList );
    
    inline void clearReferences()
    {
        treeSearch.clearReferences();
    }
    
    inline int numberOfFaces() const
    {
        return nelements;
    }
    
    inline void clearAll() 
    {
        treeSearch.clearAll();
	
	if( faces )
	{

	    for( int i = 0; i < npatches(); ++i )
	    {
	       delete[] faces[i];
	       delete[] demFaceAreas[i];
	    } 
	       
	    delete[] faces;
	    npatch_ = 0;
	    
	    delete[] nfaces_;
	    delete[] demFaceAreas;
	    nfaces_ = NULL;
	    faces = NULL;
	    demFaceAreas = 0;
	}
	
    }
           
};

#endif


