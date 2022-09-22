
#include <vector>
#include "fvCFD.H"
#include <cstdio>
#include "boundaryElementFinder.H"

// -- Face Class --

boundaryElementFinder::Face::Face() :
patchI( 0 ),
faceI( 0 )
{}


boundaryElementFinder::Face::Face( Foam::label patchI_, Foam::label faceI_ ) :
patchI( patchI_ ),
faceI( faceI_ )
{}

boundaryElementFinder::Face::~Face()
{}

void boundaryElementFinder::Face::setIDs( Foam::label patchI_, Foam::label faceI_ )
{
   patchI = patchI_;
   faceI = faceI_;
}


bool boundaryElementFinder::Face::isInside( const Foam::vector& point, bool debug ) const
{       
    // -- face is incomplete --
    if( points.size() < 3 ) return false;
	
    // -- face must have at least three points --
    Foam::vector p1 = points.at(0);
    Foam::vector p2 = points.at(1);
    Foam::vector p3 = points.at(2);

    // -- spanning vectors of the plane where face belongs --
    
    Foam::vector e1 = (p2-p1)/mag(p2-p1);
    Foam::vector e2 = (p3-p1)/mag(p3-p1);
    
    Foam::vector e3 = e1 ^ e2;

    e3 /= Foam::mag( e3 );

    // -- distance from the place --

    Foam::scalar dist = e3 & (point-p1);
    
    Foam::scalar distExtra = 0;
    
    // -- if there are more points use the extra points to determine the face accuracy --
    for( int ii = 3; ii < points.size(); ++ii )
    {
	Foam::vector cc = (points.at(ii)-p1);
	Foam::scalar distExtra_ = abs( e3 & cc );
	distExtra = distExtra_ > distExtra ? distExtra_ : distExtra;
    }
    
    
    if( debug )
    {
        printf( "dist = %e    distExtra = %e \n ", dist, distExtra );
    }
    
    // -- point is not in the place --
    if( abs( dist ) > distanceThreshold + distExtra ) return false;


    // -- check if the point is on the face --

    // -- create orthogonal basis for the plane --
    Foam::vector ee2 = e3 ^ e1;
    
    // -- matrix presenting the orthogonal basis for the face --
    Foam::tensor base( 
	     	       e1.x(), e1.y(), e1.z(), 
	     	       ee2.x(), ee2.y(), ee2.z(), 
		       e3.x(), e3.y(), e3.z() 
		     );

    std::vector<Foam::vector> projectedFace;

    // -- map the face and the point of interest on the plane --
    Foam::vector projectedPoint = base & (point-p1);
    projectedPoint.z() = 0;
            
    for( int i = 0; i < points.size(); ++ i )
    {
    	Foam::vector cc = base & ( points.at(i) - p1 );
	
	// -- remove numerical remnants from the projected points (should lie in the plane) --
	cc.z() = 0;
	projectedFace.push_back( cc );
	
    }

    Foam::scalar angleSum = 0;

    // -- start traversing the face from the last point --
    Foam::vector pp = projectedFace.at( points.size() - 1 ) - projectedPoint;

    for( int i = 0; i < points.size(); ++i )
    {
	Foam::vector pn = projectedFace.at(i) - projectedPoint;	         
	Foam::scalar vangle = vectorAngle( pp, pn );
	
	if( debug ) printf( "Angle %d  =  %f \n ", i, 360 * vangle/(fullAngle) );
	
	// -- add angle --
	angleSum += vangle;
	pp = pn;
	
    }
        
    if( debug ) printf( "Total Angle  =  %f \n ", 360 * angleSum/(fullAngle) );
    
    // -- if angles sum up to full round point is in the face --
    if( 
    	Foam::mag( angleSum - fullAngle ) < angleThreshold || 
	Foam::mag( angleSum + fullAngle ) < angleThreshold
      )
    {
	return true;
    }else{
	return false;
    }

}


Foam::scalar boundaryElementFinder::Face::vectorAngle( const Foam::vector& a, const Foam::vector& b ) const
{
    
    Foam:vector out = ( a ^ b );
    Foam::scalar in = ( a & b )/ ( Foam::mag(a) * Foam::mag(b) );
    
    if( out.z() > 0 ) // -- counter clockwise --
    {
        return Foam::acos( in );
    }else{ // -- clockwise --
        return - Foam::acos( in );
    }

}

void boundaryElementFinder::Face::addReference( int imesh, int ielement )
{
    meshList.push_back( imesh );
    elementList.push_back( ielement );    
}
	
void boundaryElementFinder::Face::getBounds( double* min, double* max ) const
{

    min[0] = INF;
    min[1] = INF;
    min[2] = INF;

    max[0] = -INF;
    max[1] = -INF;
    max[2] = -INF;	   

    for( int i = 0; i < points.size(); ++i )
    {

	const Foam::vector& p = points.at(i);

	for( int ii = 0; ii < 3; ++ii )
	{
	    if( p[ii] < min[ii] ) min[ii] = p[ii];
	    if( p[ii] > max[ii] ) max[ii] = p[ii];
	}

    }

}

void boundaryElementFinder::Face::clearReferences()
{
    meshList.clear();
    elementList.clear();
}

// -- Container Class --

boundaryElementFinder::Container::Container( double* lower_, double* upper_ )
{
    this->set( lower_, upper_ );

    // -- set tree elements to NULL pointers --
    for( int i = 0; i < 8; ++i )
       octoTree[i] = NULL;

}

boundaryElementFinder::Container::~Container()
{

    for( int i = 0; i < 8; ++i )
       if( octoTree[i] ) delete octoTree[i];

}

bool boundaryElementFinder::Container::findFace( const Foam::vector& point, Foam::label& faceId, Foam::label& patchId )
{

    int index = pointIndex( point );


    // -- check if point is on the faces in the container --
    for( int i = 0; i < faces.size(); ++i )
    {
	boundaryElementFinder::Face* ff = faces.at(i);

	if( ff->isInside( point ) ) // -- face found -- 
	{		     
	    faceId = ff->getFaceID();
	    patchId = ff->getPatchID();   
	    return true;
	}
    }

    if( index < 0 ) return false;

    boundaryElementFinder::Container* subcontainer = octoTree[index];

    if( subcontainer )
    {
	return subcontainer->findFace( point, faceId, patchId );
    }	     
    else return false;

}

bool boundaryElementFinder::Container::addReference( const Foam::vector& point, int imesh, int ielement )
{
    int index = pointIndex( point );

    // -- check if point is on the faces in the container --

    for( int i = 0; i < faces.size(); ++i )
    {
	boundaryElementFinder::Face* ff = faces.at(i);

	if( ff->isInside( point ) ) // -- face found -- 
	{	
	    ff->addReference( imesh, ielement );
	    return true;
	}
    }

    if( index < 0 ) return false;

    boundaryElementFinder::Container* subcontainer = octoTree[index];

    if( subcontainer )
    {
	return subcontainer->addReference( point, imesh, ielement );
    }else{
        return false;
    }

}

// -- returns list of all the DEM side mesh elements that belong to the fluid mesh element containing the "point" --
bool boundaryElementFinder::Container::getDEMreferences( const Foam::vector& point, std::vector<int>*& meshList, std::vector<int>*& elementList )
{

    int index = pointIndex( point );

    // -- check if point is on the faces in the container --

    for( int i = 0; i < faces.size(); ++i )
    {
	boundaryElementFinder::Face* ff = faces.at(i);

	if( ff->isInside( point ) ) // -- face found -- 
	{		     
	    meshList = ff->getMeshList();
	    elementList = ff->getElementList();
	    return true;
	}
    }

    if( index < 0 ) return false;

    boundaryElementFinder::Container* subcontainer = octoTree[index];

    if( subcontainer )
    {
	return subcontainer->getDEMreferences( point, meshList, elementList );
    }else{
	return false;
    }

}

void boundaryElementFinder::Container::getFace( const Foam::vector& point, Foam::label patchI, Foam::label faceI )
{

    int index = pointIndex( point );


    // -- check if point is on the faces in the container --
    for( int i = 0; i < faces.size(); ++i )
    {
	boundaryElementFinder::Face* ff = faces.at(i);
	
	if( patchI == ff->getPatchID() &&  faceI == ff->getFaceID() ) // -- face found -- 
	{		     
	    ff->isInside( point, true );  
	    return;
	}
    }

    if( index < 0 ) return;

    boundaryElementFinder::Container* subcontainer = octoTree[index];

    if( subcontainer )
    {
	subcontainer->getFace( point, patchI, faceI );
    }	     
    else return;

}

// -- add face in the containers -- return fase if face does not fit in the container
bool boundaryElementFinder::Container::addFace( boundaryElementFinder::Face& ff )
{

    // -- check if face is contained in the container --
    if( !isContained( ff ) ) return false;

    this->addNewFace( ff );

    return true;

} 

void boundaryElementFinder::Container::clearAll()
{
    faces.clear();

    for( int i = 0; i < 8; ++i )
       if( octoTree[i] ) octoTree[i]->clearAll();

}

void boundaryElementFinder::Container::clearReferences()
{ 

    for( int i = 0; i < faces.size(); ++i )
    {
	faces.at(i)->clearReferences();	         
    }

    for( int i = 0; i < 8; ++i )
       if( octoTree[i] ) octoTree[i]->clearReferences();

}

void boundaryElementFinder::Container::addNewFace( boundaryElementFinder::Face& ff )
{

    // -- check if face is contained in any of the sub-containers --
    int index = containedIndex( ff );

    if( index >= 0 )
    {

	Container* subcontainer = octoTree[index];

	if( subcontainer ) subcontainer->addNewFace( ff );
	else
	{
	   // -- create new sub-container --
	   int ii[3] = {0,0,0};
	   int index_ = index;


	   ii[0] = index%2;

	   index = ( index - ii[0] )/2;
	   ii[1] = index%2;
	   index = ( index - ii[1] )/2;

	   ii[2] = index%2;

	   double newLower[3];
	   double newUpper[3];

	   for( int i = 0; i < 3; ++i )
	   {
	       newLower[i] = ii[i] == 0 ? lower[i] : mid[i];
	       newUpper[i] = ii[i] == 0 ? mid[i] : upper[i];
	   }

	   octoTree[index_] = new Container( newLower, newUpper );
	   octoTree[index_]->addNewFace( ff );

	}

    }else // -- face is not contained in the sub-containers -> add face in the present container
    {
	faces.push_back( &ff );
    }

}

void boundaryElementFinder::Container::set( double* lower_, double* upper_ )
{
    for( int ii = 0; ii < 3; ++ii )
    {
       lower[ii] = lower_[ii];
       upper[ii] = upper_[ii];
       mid[ii] = (upper[ii] + lower[ii])/2.0; 
    }
}

// -- check if the face is fully enclosed in the container --
bool boundaryElementFinder::Container::isContained( const Face& face ) const
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

// -- return index of the octoTree member, which contains face -- returns -1 if not contained in any of the sub-containers
int boundaryElementFinder::Container::containedIndex( const Face& face ) const
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

// -- return index of the sub-container, -1 if point is not in the subcontainer
int boundaryElementFinder::Container::pointIndex( const Foam::vector& point ) const
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

int boundaryElementFinder::Container::countEmptyFaces() const
{
   int nn = 0;
   
   for( int i = 0; i < faces.size(); ++i )
   {
       if( faces.at(i)->nreferences() == 0 )
          ++nn;
   }

   for( int i = 0; i < 8; ++i )
       if( octoTree[i] ) nn += octoTree[i]->countEmptyFaces();   
   
   return nn;
   
}

/////////////////////////// -- Boundary Element Finder class -- //////////////////////

boundaryElementFinder::boundaryElementFinder( double* lower_, double* upper_ ) :
nelements( 0 ),
treeSearch( lower_, upper_ ),
faces( NULL ),
demFaceAreas( NULL ),
npatch_( 0 ),
nfaces_( 0 )
{}

boundaryElementFinder::~boundaryElementFinder()
{
    this->clearAll();   
}

// OpenFoam faces
bool boundaryElementFinder::addFace( const Foam::face& face_, const Foam::pointField& face_points, Foam::label patchI, Foam::label faceI )
{

    boundaryElementFinder::Face& ff = faces[patchI][faceI];
    
    ff.setIDs( patchI, faceI );
    
    forAll( face_, facei )
    {
	Foam::vector vertice = face_points[facei];
	ff.addPoint( vertice );
    }

   bool faceAdded = treeSearch.addFace( ff );
   if( faceAdded ) ++nelements;
   return faceAdded;

}

void boundaryElementFinder::initPatches( int npatch )
{
    this->npatch_ = npatch;
    nfaces_ = new int[npatches()];
    demFaceAreas = new double*[npatches()];
    
    for( int i = 0; i < npatches(); ++i )
       nfaces_[i] = 0;
       
    faces = new boundaryElementFinder::Face*[npatches()];
    
}

void boundaryElementFinder::initPatch( int iPatch, int nface )
{
    faces[iPatch] = new boundaryElementFinder::Face[nface];
    nfaces_[iPatch] = nface;
    demFaceAreas[iPatch] = new double[nface];
    
    for( int i = 0; i < nface; ++i )
        demFaceAreas[iPatch][i] = 0;
} 

bool boundaryElementFinder::getDEMreferences( Foam::label patchI, Foam::label faceI, std::vector<int>*& meshList, std::vector<int>*& elementList )
{
    boundaryElementFinder::Face& ff = faces[patchI][faceI];
    meshList = ff.getMeshList();
    elementList = ff.getElementList();   
    
    if( elementList->size() > 0 && meshList->size() > 0 ) 	return true;
    else							return false;
    
}


