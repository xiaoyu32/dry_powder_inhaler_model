/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.


    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "twoWayMPIpu.H"
#include "addToRunTimeSelectionTable.H"
#include "clockModel.H"
#include "pair.h"
#include "force.h"
#include "forceModel.H"
#include "face.H"

#define PARTICLE_ALLOCATION_FRACTION 0.5

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoWayMPIpu, 0);

addToRunTimeSelectionTable
(
    dataExchangeModel,
    twoWayMPIpu,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



// Construct from components
twoWayMPIpu::twoWayMPIpu
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dataExchangeModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    bufferSizeForAlloc_(0)
{
    // set max nr of particles from dict
    //Info << "twoWayMPIpu.C- this should no longer be needed" << endl;
    //maxNumberOfParticles_ = readScalar(propsDict_.lookup("maxNumberOfParticles"));

    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    if (me < nprocs) liggghts = 1;
    else liggghts = MPI_UNDEFINED;

    MPI_Comm_split(MPI_COMM_WORLD,liggghts,0,&comm_liggghts);

    // Sub-domain min&max coordinates
    // const pointField& pp = sm.mesh().points();
    
    vector imin(1.e10,1.e10,1.e10);
    vector imax(-1.e10,-1.e10,-1.e10);
    vector unit(0,0,0);
    
    forAll(sm.mesh().boundaryMesh(), patchI)
    {

       const List<Foam::face> faces = sm.mesh().boundaryMesh()[patchI].localFaces();    
       pointField meshPoints = sm.mesh().boundaryMesh()[patchI].localPoints(); 
       pointField face_normals = sm.mesh().boundaryMesh()[patchI].faceNormals();

       forAll( faces, ii )
       {

	   Foam::face face_ = faces[ii];
	   
	   pointField face_points = face_.points( meshPoints );    

	   forAll( face_, facei )
	   { 
	       vector pC = face_points[facei];
	       
	       for(int i=0; i<3; i++)
               {

		  if(i==0) unit=vector(1,0,0);
        	  if(i==1) unit=vector(0,1,0);
        	  if(i==2) unit=vector(0,0,1);

        	  scalar iminUpdate = (pC & unit);
        	  scalar imaxUpdate = (pC & unit);

		  // -- make the domains slightly larger to avoid particles not sent to any processor --
        	  if(iminUpdate<imin[i]) imin[i] = iminUpdate;
        	  if(imaxUpdate>imax[i]) imax[i] = imaxUpdate;

               }
	       
	   }

       }

      //const polyPatch& patchFound = sm.mesh().boundaryMesh()[patchI];
      //labelList labelPatchFound(patchFound.meshPoints() );
      //pointField meshPoints(sm.mesh().points());
      
    }


    // Min, max x-coordinates
    //double coordCPUs[nprocs*6]; 
    double coordCPUs[6];
    double allCoordCPUs[nprocs*6]; 

    //if ( me == 0) { 
       //allCoordCPUs = (double *)malloc(nprocs*6*sizeof(double)); 
    //} 

    // Min, max x-coordinates       
    coordCPUs[0] = imin[0] - asciiPrecission; //min(pp & vector(1,0,0));
    coordCPUs[1] = imax[0] + asciiPrecission; //max(pp & vector(1,0,0));

    // Min, max y-coordinates               
    coordCPUs[2] = imin[1] - asciiPrecission; //min(pp & vector(0,1,0));
    coordCPUs[3] = imax[1] + asciiPrecission; //max(pp & vector(0,1,0));

    // Min, max z-coordinates               
    coordCPUs[4] = imin[2] - asciiPrecission; //min(pp & vector(0,0,1));
    coordCPUs[5] = imax[2] + asciiPrecission; //max(pp & vector(0,0,1));

    MPI_Gather(coordCPUs, 6, MPI_DOUBLE, allCoordCPUs, 6, MPI_DOUBLE, 0,MPI_COMM_WORLD);  

    if (me == 0) 
    {

      //const fileName cpusPath(propsDict_.lookup("cpusPath"));
      const fileName cpusPath("../DEM/in.cpus");
      char * cpusPathChar = (char*)cpusPath.c_str();

      Info<<"Creating cpus coords file '"<< cpusPath.c_str() <<"'"<<endl;


      //set file pointer
      ofstream outputPtr(cpusPathChar);


      if (!outputPtr) {
        printf("ERROR: Could not open cpu coords file \n");
        MPI_Abort(MPI_COMM_WORLD,1);
      }
      // -- use scientific syntax --
      outputPtr << std::scientific;
      

      outputPtr << nprocs << nl;
      for(int index = 0; index < nprocs*6; ++index)
      {
        outputPtr << allCoordCPUs[index] << nl; 
      }

}
    
    write_boundary_info( sm );
    

    Info<<"Starting up LIGGGHTS for first time execution"<<endl;

    // open LIGGGHTS input script
    FILE *fp=NULL;
    if (me == 0)
    {
      // read path from dictionary
      const fileName liggghtsPath(propsDict_.lookup("liggghtsPath"));
      char * liggghtsPathChar = (char*)liggghtsPath.c_str();

      Info<<"Executing input script '"<< liggghtsPath.c_str() <<"'"<<endl;

      fp = fopen(liggghtsPathChar,"r");

      if (fp == NULL) {
        printf("ERROR: Could not open LIGGGHTS input script\n");
        MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    
    if (liggghts == 1) lmp = new LAMMPS_NS::LAMMPS(0,NULL,comm_liggghts);

    int n;
    char line[1024];
    while (1) {
      if (me == 0) {
        if (fgets(line,1024,fp) == NULL) n = 0;
        else n = strlen(line) + 1;
        if (n == 0) fclose(fp);
      }
      MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
      if (n == 0) break;
      MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
      if (liggghts == 1) lmp->input->one(line);
    }

    // get DEM time step size
    DEMts_ = lmp->update->dt;
    Info<< "DEMts_ = lmp->update->dt "<< DEMts_<<endl;
    checkTSsize();
			
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoWayMPIpu::~twoWayMPIpu()
{}

void twoWayMPIpu::write_boundary_info( cfdemCloud& sm )
{

    int me, nprocs;
    int n;
    int np;
    
    char * filename = "../DEM/in.cpus.faces";
    double buffer[50];
    
    MPI_Status status;
    MPI_Comm_rank( MPI_COMM_WORLD, &me );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    
    if( me == 0 )
    {
        //write the root boundary face informations
	ofstream file;
	file.open( filename );
	
	// -- use scientific syntax --
	file<< std::scientific;
	
	file<<"Processor: "<<me<<nl;
	
	forAll(sm.mesh().boundaryMesh(), patchI)
        {

            //const polyPatch& patchFound = sm.mesh().boundaryMesh()[patchI];
            //labelList labelPatchFound(patchFound.meshPoints() );
            //pointField meshPoints(sm.mesh().points());

	    const List<Foam::face> faces = sm.mesh().boundaryMesh()[patchI].localFaces();    
	    pointField meshPoints = sm.mesh().boundaryMesh()[patchI].localPoints(); 
	    pointField face_normals = sm.mesh().boundaryMesh()[patchI].faceNormals();
	    
	    forAll( faces, ii )
            {
		
		Foam::face face_ = faces[ii];
		vector normal = face_normals[ii];
		
		double r = normal[0]*normal[0] +
		           normal[1]*normal[1] +
			   normal[2]*normal[2];
		
		
	        file<< "Normal: " << normal[0] 
		    << " , " << normal[1] 
		    << " , " << normal[2] <<nl;
		    
		pointField face_points = face_.points( meshPoints );    
		    
		forAll( face_, facei )
		{ 
		    vector vertice = face_points[facei];
		    file<< vertice[0] << " , "<< vertice[1] << " , "<< vertice[2] << "  ";
		}
		
		file<<" "<<nl;
		
	    }
	
	}
	
	//Receive boundary face information from neighboring processors
	for( int id = 1; id < nprocs; ++id )
	{
	    //number of boundary faces
	    MPI_Recv( &n, 1, MPI_INT, id, 99, MPI_COMM_WORLD, &status );
 
	    file<<"Processor: "<<id<<nl;
	    
	    for( int i = 0; i < n; ++i )
	    {
	       
	        MPI_Recv( buffer, 3, MPI_DOUBLE, id, 99, MPI_COMM_WORLD, &status );
//                       MPI_Barrier( MPI_COMM_WORLD );
		

		file<< "Normal: " << buffer[0] 
		    << " , " << buffer[1] 
		    << " , " << buffer[2] <<nl;
	       
	        MPI_Recv( &np, 1, MPI_INT, id, 99, MPI_COMM_WORLD, &status );

		MPI_Recv( buffer, 3*np, MPI_DOUBLE, id, 99, MPI_COMM_WORLD, &status );
//		        MPI_Barrier( MPI_COMM_WORLD );


	        for( int j = 0; j < np; ++j )
		{
		    file<< buffer[3*j] << " , " << buffer[3*j+1] << " , " << buffer[3*j+2] << "  ";    
		}
		
		file<<" "<<nl;
	       
	    }
	    
	    
	}
        

        file.close();


	
    }else{
        
	std::vector< vector > ninfo;
	std::vector< std::vector< vector > > finfo;
	std::vector< vector > buff;
	
	//construct vectors consisting of vertices and face normals
	forAll(sm.mesh().boundaryMesh(), patchI)
        {

            //const polyPatch& patchFound = sm.mesh().boundaryMesh()[patchI];
            //labelList labelPatchFound(patchFound.meshPoints() );
            //pointField meshPoints(sm.mesh().points());

	    const List<Foam::face> faces = sm.mesh().boundaryMesh()[patchI].localFaces();    
	    pointField meshPoints = sm.mesh().boundaryMesh()[patchI].localPoints();
	    pointField face_normals = sm.mesh().boundaryMesh()[patchI].faceNormals();
	    
	    forAll( faces, ii )
            {
		
		Foam::face face_ = faces[ii];
		vector normal = face_normals[ii];
	        ninfo.push_back( normal );
		    
		pointField face_points = face_.points( meshPoints );    
		
		buff.clear();
		    
		forAll( face_, facei )
		{ 
		    vector vertice = face_points[facei];
		    buff.push_back( vertice );
		}
		
		finfo.push_back( buff );
	    }
	
	}
	
	
	//Send data to root
	n = (int)ninfo.size();
	 
	MPI_Send( &n, 1, MPI_INT, 0, 99, MPI_COMM_WORLD ); 
 
	for( int i = 0; i < ninfo.size(); ++i )
	{
	    
	    for( int j = 0; j < 3; ++j )
	       buffer[j] = ninfo.at(i)[j];
	    
	    MPI_Send( buffer, 3, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD );   
        //            MPI_Barrier( MPI_COMM_WORLD );

	    np = (int)finfo.at(i).size();
	    
	    MPI_Send( &np, 1, MPI_INT, 0, 99, MPI_COMM_WORLD );

	    for( int j = 0; j < np; ++j )
	       for( int k = 0; k < 3; ++k )
	       {
	           buffer[3*j+k] = finfo.at(i).at(j)[k];
	       }
	    
	    MPI_Send( buffer, 3*np, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD );   
          //               MPI_Barrier( MPI_COMM_WORLD );

	}
	
	 
    }
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// mostafa replaced this function

/*

char* twoWayMPIpu::wordToChar(word& inWord) const

  {
    string HH = string(inWord);
    return const_cast<char*>(HH.c_str());
}

*/


char* twoWayMPIpu::wordToChar(word& inWord) const
{
    return const_cast<char*>(inWord.c_str());
}



// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
void twoWayMPIpu::getData
(
    word name,
    word type,
    double ** const& field,
    label step
) const
{
    char* charName = wordToChar(name);
    char* charType = wordToChar(type);
    data_liggghts_to_of(charName,charType, lmp, (void*&) field, (char *)"double");

}

void twoWayMPIpu::getData
(
    word name,
    word type,
    int ** const& field,
    label step
) const
{
    char* charName = wordToChar(name);
    char* charType = wordToChar(type);
    data_liggghts_to_of(charName,charType, lmp, (void*&) field, (char *)"int");

}

void twoWayMPIpu::giveData
(
    word name,
    word type,
    double ** const& field,
    const char* datatype
) const
{
    char* charName = wordToChar(name);
    char* charType = wordToChar(type);
    char* charDatatype= const_cast<char*> (datatype);
    data_of_to_liggghts(charName,charType, lmp, (void*) field,charDatatype);
}
//============
// double **
void Foam::twoWayMPIpu::allocateArray
(
    double**& array,
    double initVal,
    int width,
    int length
) const
{
    //if(length==-1) then LIGGGHTS uses own length data
    allocate_external_double(array, width,length,initVal,lmp);
}



void Foam::twoWayMPIpu::allocateArray
(
    double**& array,
    double initVal,
    int width,
    const char* length
) const
{
    //if(length==-1) then LIGGGHTS uses own length data
    char* charLength= const_cast<char*> (length);
    allocate_external_double(array, width,charLength,initVal,lmp);
}
void Foam::twoWayMPIpu::destroy(double** array,int len) const
{
    if (array == NULL) return;

    //for ( int i = 0; i < len; i++ ) // does not work
    for ( int i = 0; i < 1; i++ )
        free(array[i]); 

    free(array);
}
//============
// int **
void Foam::twoWayMPIpu::allocateArray
(
    int**& array,
    int initVal,
    int width,
    int length
) const
{

    //if(length==-1) then LIGGGHTS uses own length data
    allocate_external_int(array, width,length,initVal,lmp);
}

void Foam::twoWayMPIpu::allocateArray
(
    int**& array,
    int initVal,
    int width,
    const char* length
) const
{
    //if(length==-1) then LIGGGHTS uses own length data
    char* charLength= const_cast<char*> (length);
    allocate_external_int(array, width,charLength,initVal,lmp);
}
void Foam::twoWayMPIpu::destroy(int** array,int len) const
{
    if (array == NULL) return;

    //for ( int i = 0; i < len; i++ ) // does not work
    for ( int i = 0; i < 1; i++ )
        free(array[i]); 


    free(array);
}
//============
// int *
void Foam::twoWayMPIpu::destroy(int* array) const
{
    if (array == NULL) return;
    free(array);
}
//============
// double *
void Foam::twoWayMPIpu::destroy(double* array) const
{
    if (array == NULL) return;
    free(array);
    
}
//============

bool Foam::twoWayMPIpu::couple() const
{
    bool coupleNow = false;
    
    if ( doCoupleNow() )
    {
        this->coupleEnd_ = true;
        couplingStep_++;
        coupleNow = true;

        // start liggghts
        if (liggghts == 1)
        {
            /*// hardcoded run commands
            char lammpsRunCommand[80];
            if (couplingStep()==1) sprintf(lammpsRunCommand,"run %d",int(couplingInterval_));
            else           sprintf(lammpsRunCommand,"run %d pre no",int(couplingInterval_));
            Info << "old script would Executing command: '"<<lammpsRunCommand <<"'"<< endl;
            lmp->input->one(lammpsRunCommand);*/

            // run commands from liggghtsCommands dict
            Info<<"Starting up LIGGGHTS" << endl;
            particleCloud_.clockM().start(3,"LIGGGHTS");
            forAll(particleCloud_.liggghtsCommandModelList(),i)
            {
                if(particleCloud_.liggghtsCommand()[i]().runCommand(couplingStep()))
                {
                    const char* command = particleCloud_.liggghtsCommand()[i]().command();
                    Info << "Executing command: '"<< command <<"'"<< endl;
                    lmp->input->one(command);
                }
            }
            particleCloud_.clockM().stop("LIGGGHTS");
            Info<<"LIGGGHTS finished"<<endl;
        }

        // give nr of particles to cloud
        //double newNpart = liggghts_get_maxtag(lmp);
	//MPI_pu development
        int newNpart = liggghts_get_localtag(lmp, true);		
	
        setNumberOfParticles(newNpart);

        // re-allocate arrays of cloud
        particleCloud_.clockM().start(4,"LIGGGHTS_reallocArrays");
	
	if( newNpart > bufferSizeForAlloc_ || 
            newNpart < PARTICLE_ALLOCATION_FRACTION * bufferSizeForAlloc_ )
        {
	    bufferSizeForAlloc_ = newNpart;
	    particleCloud_.reAllocArrays();
	    Info << "in CFD, reallocated ..." << endl; 
	}
	
        particleCloud_.clockM().stop("LIGGGHTS_reallocArrays");  			
	
    }

    return coupleNow;
}

int Foam::twoWayMPIpu::getNumberOfParticles() const
{
    return liggghts_get_localtag(lmp, true);
}

int Foam::twoWayMPIpu::getNumberOfClumps() const
{
    Warning << "liggghts_get_maxtag_ms(lmp) is commented here!" << endl;
    return -1;
//    return liggghts_get_maxtag_ms(lmp);
}

void* Foam::twoWayMPIpu::getLmp() const
{
   return (void*)lmp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
