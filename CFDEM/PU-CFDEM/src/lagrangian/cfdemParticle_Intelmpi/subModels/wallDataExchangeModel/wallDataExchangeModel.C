#include "mpi.h"
#include "fvCFD.H"
#include "cfdemCloud.H"
#include "wallDataExchangeModel.H"
#include "boundaryElementFinder.H" 

#include <cstdio>
#include <library_cfd_coupling.h>



using namespace Foam;

wallDataExchangeModel::wallDataExchangeModel( void *lmp_, cfdemCloud& sm ) :
lmp( lmp_ ),
particleCloud_( sm ),
nodeCenters( "center" ),
faces( "node" ),
finderCFD( NULL )
{
   
   // -- construct boundary element finder --
   
   
   double min[3] = {INF,INF,INF};
   double max[3] = {-INF,-INF,-INF}; 
   
   
   forAll(particleCloud_.mesh().boundaryMesh(), patchI)
   {
      const polyPatch& patchFound = particleCloud_.mesh().boundaryMesh()[patchI];
      labelList labelPatchFound(patchFound.meshPoints() );
      pointField meshPoints(particleCloud_.mesh().points());
      
      forAll(labelPatchFound , label)
      {
         vector pC = meshPoints[labelPatchFound[label]];
         
	 for(int i=0; i<3; i++)
         {
            
	    if( pC[i] < min[i] ) min[i] = double( pC[i] ) - 1e-8; 
	    if( pC[i] > max[i] ) max[i] = double( pC[i] ) + 1e-8;
	  
         }
	 
      }
         
   }
   
   // -- initialize boundary element finder --
   finderCFD = new boundaryElementFinder( min, max );
   finderCFD->initPatches( particleCloud_.mesh().boundaryMesh().size() );
   
   forAll(particleCloud_.mesh().boundaryMesh(), patchI)
   {
       finderCFD->initPatch( patchI, particleCloud_.mesh().boundaryMesh()[patchI].localFaces().size() );       
   }
   
   
   // -- insert boundary faces into the finder --
   forAll(particleCloud_.mesh().boundaryMesh(), patchI)
   {
       const List<Foam::face> faces = particleCloud_.mesh().boundaryMesh()[patchI].localFaces();
       pointField meshPoints = particleCloud_.mesh().boundaryMesh()[patchI].localPoints();
       
       forAll( faces, ii )
       {
           
	   Foam::face face_ = faces[ii];
	   pointField face_points = face_.points( meshPoints );
	   
	   // -- add face to the finder -- 
	   finderCFD->addFace( face_, face_points, patchI, ii );
	   	   
       }
   
   }
   
   bool& wallUpdate = isWallUpdated();
   wallUpdate = true;
   
   this->update();
   
}



wallDataExchangeModel::~wallDataExchangeModel()
{
   delete finderCFD;
}

void wallDataExchangeModel::update()
{
      
   bool& wallUpdate = isWallUpdated();
   
   // -- DEM wall is upto date --
   if( !wallUpdate ) return; 
   
   Info<<"Wall exchange model updated."<<endl;
   
   
   wallUpdate = false; // -- mark wall exchange model updated --
      
   // -- get DEM node centers --
   getData( nodeCenters );
   
   // -- get DEM face vertices --
   getData( faces );
   
   finderCFD->clearReferences();
      
   nodeCenters.begin();
   
   double* pos = NULL;
   int imesh;
   int ielement;
   
   int refCounter = 0;
   
   // -- insert indices of corresponding DEM mesh elements to CFD mesh elements --
   do
   { 
	pos = nodeCenters.nextVector(imesh,ielement);
	if( !pos ) break;
	
	Foam::vector location( pos[0],pos[1],pos[2] );	
	
	if( finderCFD->addReference( location, imesh, ielement ) )
	   ++refCounter;	
       
   }while( pos );
   
   
   if( refCounter < nodeCenters.size() )
   {
       printf( "Warning: No matching CFD boundary element found for %d/%d DEM wall elements! \n", nodeCenters.size()-refCounter, nodeCenters.size() );
   }
   
   // -- calculate total DEM face areas for each CFD face --
   double**** faceData = faces.ref();
   
   
   forAll(particleCloud_.mesh().boundaryMesh(), patchI)
   {
       const pointField& face_centers = particleCloud_.mesh().boundaryMesh()[patchI].faceCentres();
       
       forAll( face_centers, ii )
       {
           
	   finderCFD->faceArea( patchI, ii ) = 0;
	   
	   std::vector<int>* meshList = NULL;
	   std::vector<int>* elementList = NULL;
	   
	   if( !finderCFD->getDEMreferences( patchI, ii, meshList, elementList ) )
	   {
	       // -- skip faces that do not have DEM counterpart --
               continue;
	   }
	   
	   for( int i = 0; i < meshList->size(); ++i )
	   {
	        Foam::scalar demFaceArea = computeArea( (const double**)faceData[ meshList->at(i) ][elementList->at(i)], faces.getNvec() );
		finderCFD->faceArea( patchI, ii ) += demFaceArea;
	   }
	   
       }
   }
   
   /*
   int emptyReferences = finderCFD->countEmptyFaces();
   
   if( emptyReferences > 0 )
   {
       printf( "Warning: there are %d/%d faces without DEM wall element reference! \n", emptyReferences, finderCFD->numberOfFaces() );
   }
   */
   
}

void wallDataExchangeModel::initData( WallValueContainer& container )
{
   int len = 0;
   int nvec = 0; 
   int nmesh = nMeshes( lmp );
   int nelements = 0;
   
   if( nmesh <= 0 ) return; // -- there are no mesh objects on liggghts --
   
   this->update();
   
   container.allocateData( nmesh );
   
   for( int imesh = 0; imesh < nmesh; ++imesh )
   {
       nelements = numberOfElements( lmp, imesh );
              
       dataSize( lmp, container.getID(), len, nvec );
          
       //printf( "ID = %s      imesh = %d   nelements = %d    len = %d     nvec = %d  \n", container.getID(), imesh,nelements,len,nvec );	   //fixme
       //MPI_Barrier( MPI_COMM_WORLD );
       	  
       if( nelements > 0 ) container.allocateMesh( imesh, nelements, len, nvec );
       
       // -- could not locate the property (this is OK) --
       if( len == 0 && nvec == 0 )
       {
           //Info<<"Warning: Failed to locate DEM wall property '"<<container.getID()<<"' for mesh '"<<imesh<<"'. Check the DEM settings."<<endl;
       }
       
   }
   
}

void wallDataExchangeModel::getData( WallValueContainer& container )
{
   initData( container );
   wall_data_liggghts_to_of( container.getID(), container.ref(), this->lmp );
}


void wallDataExchangeModel::giveData( WallValueContainer& container )
{
    this->update();
    wall_data_of_to_liggghts( container.getID(), container.ref(), this->lmp );
}

Foam::autoPtr<wallDataExchangeModel> wallDataExchangeModel::New
(
    void* lmp,
    cfdemCloud& sm
)
{
    Info<<"Setting up wall data exchange model..."<<endl;
    return autoPtr<wallDataExchangeModel>
    	   (
      		new wallDataExchangeModel( lmp, sm )
           );
}



Foam::scalar wallDataExchangeModel::computeArea( const double** face, int nvec ) const
{
    
    if( !face ) return 0;
    
    Foam::scalar area = 0;
    
    // -- face must have atleast 3 points --
    if( nvec <= 2 ) return area;
    
    Foam::vector p1( face[0][0], face[0][1], face[0][2] ); 
    Foam::vector p2( face[1][0], face[1][1], face[1][2] );
    
    for( int i = 2; i < nvec; ++i )
    {
        
	Foam::vector p3( face[i][0], face[i][1], face[i][2] );
	
	Foam::vector pA = p2-p1;
	Foam::vector pB = p3-p1;
	
	area += Foam::mag( pA ^ pB )/2.0;
	p2 = p3;
    }
    
    return area;
    
}



/* 
 *    **********  SURFACE FIELDS ************
 */
void wallDataExchangeModel::transferDataToField( WallValueContainer& container, Foam::surfaceScalarField& field )
{
    
    Info<<"Transferring data from DEM wall container to surfaceScalarField..."<<endl;
    field.internalField() = 0;    
    
    forAll( field.boundaryField(), patchI )
    {
	        
	forAll( field.boundaryField()[patchI], faceI )
	{
	    field.boundaryField()[patchI][faceI] = 0;
	}
	
    }

    double* pos = NULL;
    
    // -- begin iterators --
    nodeCenters.begin();
    container.begin();
    faces.begin();
      
    int failedCounter = 0;
    int faceCounter = 0;  
        
	
    do
    {
       
       double value = container.nextScalar();
              
       pos = nodeCenters.nextVector();
       
       Foam::scalar demFaceArea = computeArea( (const double**)faces.nextElement(), faces.getNvec() );
	      
       if( !pos ) break;
       //printf( "value = %e   counter = %d \n", value, counter++ );
       
       ++faceCounter;
       
       Foam::vector location( pos[0], pos[1], pos[2] );
       
       Foam::label faceI = -1;
       Foam::label patchI = -1;
       
       // -- find corresponding face ID and patch 
       if( !finderCFD->findFace( location, faceI, patchI ) )
       {
           ++failedCounter;
           continue;
       }
       
       Foam::scalar cfdFaceArea = finderCFD->faceArea( patchI, faceI ); // -- total sum of DEM face areas (ensures conservation) --
       Foam::scalar& cfdFaceValue = field.boundaryField()[patchI][faceI];
       
       cfdFaceValue += value * ( demFaceArea / cfdFaceArea ); 
       //cfdFaceValue = value;      
	     
       //printf( "cfdFaceValue = %e    demFaceArea = %e    cfdFaceArea = %e \n", cfdFaceValue, demFaceArea, cfdFaceArea );
       
       
    }while( pos );
    
    if( failedCounter > 0 ) printf( "Warning: %d/%d faces not found! \n", failedCounter, faceCounter );
    
    
    //int counter = 0;
    //int setCounter = 0;
    
    /*
    forAll( field.boundaryField(), patchI )
    {
	        
	forAll( field.boundaryField()[patchI], faceI )
	{
	    ++counter;
	    if( field.boundaryField()[patchI][faceI] < 0 )
	    {
	        ++setCounter;
	    }
	}
	
    }
    
    if( setCounter > 0 ) printf( "Warning: %d/%d face values not set! \n", setCounter, counter );
    */
}




void wallDataExchangeModel::transferDataFromField( WallValueContainer& container, Foam::surfaceScalarField& field )
{
    
   double**** containerData = container.ref(); 
    
   int counter = 0; 
    
   forAll(field.mesh().boundaryMesh(), patchI)
   {
       const pointField& face_centers = field.mesh().boundaryMesh()[patchI].faceCentres();
       
       forAll( face_centers, ii )
       {
           
	   std::vector<int>* meshList = NULL;
	   std::vector<int>* elementList = NULL;
	   
	   if( !finderCFD->getDEMreferences( patchI, ii, meshList, elementList ) )
	   {
	       // -- the CFD face does not have DEM counter part --
	       continue;
	   }
	   
	   if( !meshList ){ // -- this should never happen --
	      printf( "Error: retrieved NULL reference (wallExchangeModel.C line 358)!\n" );
	      continue;
	   }
	   
	   Foam::scalar& faceValue = field.boundaryField()[patchI][ii];
	   
	   // -- assign each DEM mesh cell to the same surface density as the CFD cell --
	   for( int i = 0; i < meshList->size(); ++i ){
	         ++counter;
	         containerData[meshList->at(i)][elementList->at(i)][0][0] = faceValue; 
	   }
       }
   
   }
   
   if( counter < container.size() ) printf( "Warning: Values not set for %d/%d faces DEM faces.\n", container.size()-counter, container.size() );
   
   
}

/* 
 *    **********  VOLUME FIELDS ************
 */
 
void wallDataExchangeModel::transferDataToField( WallValueContainer& container, Foam::volScalarField& field )
{
    
    Info<<"Transferring data from DEM wall container to volScalarField..."<<endl;
    field.internalField() = 0;    
    
    forAll( field.boundaryField(), patchI )
    {
	        
	forAll( field.boundaryField()[patchI], faceI )
	{
	    field.boundaryField()[patchI][faceI] = 0;
	}
	
    }

    double* pos = NULL;
    
    // -- begin iterators --
    nodeCenters.begin();
    container.begin();
    faces.begin();
      
    int failedCounter = 0;
    int faceCounter = 0;  
        
	
    do
    {
       
       double value = container.nextScalar();
              
       pos = nodeCenters.nextVector();
       
       Foam::scalar demFaceArea = computeArea( (const double**)faces.nextElement(), faces.getNvec() );
	      
       if( !pos ) break;
       //printf( "value = %e   counter = %d \n", value, counter++ );
       
       ++faceCounter;
       
       Foam::vector location( pos[0], pos[1], pos[2] );
       
       Foam::label faceI = -1;
       Foam::label patchI = -1;
       
       // -- find corresponding face ID and patch 
       if( !finderCFD->findFace( location, faceI, patchI ) )
       {
           ++failedCounter;
           continue;
       }
       
       Foam::scalar cfdFaceArea = finderCFD->faceArea( patchI, faceI ); // -- total sum of DEM face areas (ensures conservation) --
       Foam::scalar& cfdFaceValue = field.boundaryField()[patchI][faceI];
       
       cfdFaceValue += value * ( demFaceArea / cfdFaceArea ); 

    }while( pos );
    
    if( failedCounter > 0 ) printf( "Warning: %d/%d faces not found! \n", failedCounter, faceCounter );

}

void wallDataExchangeModel::transferDataFromField( WallValueContainer& container, Foam::volScalarField& field )
{
    
   double**** containerData = container.ref(); 
    
   int counter = 0; 
    
   forAll(field.mesh().boundaryMesh(), patchI)
   {
       const pointField& face_centers = field.mesh().boundaryMesh()[patchI].faceCentres();
       
       forAll( face_centers, ii )
       {
           
	   std::vector<int>* meshList = NULL;
	   std::vector<int>* elementList = NULL;
	   
	   if( !finderCFD->getDEMreferences( patchI, ii, meshList, elementList ) )
	   {
	       // -- the CFD face does not have DEM counter part --
	       continue;
	   }
	   
	   if( !meshList ){ // -- this should never happen --
	      printf( "Error: retrieved NULL reference (wallExchangeModel.C line 358)!\n" );
	      continue;
	   }
	   
	   Foam::scalar& faceValue = field.boundaryField()[patchI][ii];
	   
	   // -- assign each DEM mesh cell to the same surface density as the CFD cell --
	   for( int i = 0; i < meshList->size(); ++i ){
	         ++counter;
	         containerData[meshList->at(i)][elementList->at(i)][0][0] = faceValue; 
	   }
       }
   
   }
   
   if( counter < container.size() ) printf( "Warning: Values not set for %d/%d faces DEM faces.\n", container.size()-counter, container.size() );
   
   
}



