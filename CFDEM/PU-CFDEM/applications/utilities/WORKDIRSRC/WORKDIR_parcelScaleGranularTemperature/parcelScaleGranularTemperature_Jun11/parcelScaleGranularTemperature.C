/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    postCalc

Description
    Generic wrapper for calculating a quantity at each time

\*---------------------------------------------------------------------------*/

#include "createParcels.H"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	
 	void superParcelScaleTemp
	(
		cfdemCloud& sm,
		double**& parcelVelocities_,
		vector& meanParcelVel_,			
		scalar& aveSuperQparcel2_
	)
	{		
		aveSuperQparcel2_ = 0 ;		
		
		for(int index = 0; index <  sm.numberOfParticles(); index++)
	        {			
			aveSuperQparcel2_ += 1./3.* (   ( parcelVelocities_[index][0] - meanParcelVel_[0] ) * ( parcelVelocities_[index][0] - meanParcelVel_[0] ) 
			                              + ( parcelVelocities_[index][1] - meanParcelVel_[1] ) * ( parcelVelocities_[index][1] - meanParcelVel_[1] )
			                              + ( parcelVelocities_[index][2] - meanParcelVel_[2] ) * ( parcelVelocities_[index][2] - meanParcelVel_[2] ) );
		}												
	}
		

	template <typename T>
	string NumberToString ( T Number )
  	{
     		OStringStream ss; //ostringstream ss;
     		ss << Number;
     		return ss.str();
  	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{	

	Foam::timeSelector::addOptions();
	#   include "addRegionOption.H"
	Foam::argList::addBoolOption
	(
	"noWrite",
	"suppress writing results"
	);
	#include "addDictOption.H"

	#include "setRootCase.H"
	#include "createTime.H"
	Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
	#include "createNamedMesh.H"

	#include "createFields.H" 	 

	// Create particle cloud			
	cfdemCloud particleCloud(mesh);	
	particleCloud.reAllocArrays();

	double **positions;
	double **velocities;
	double **radii;
	double **cellID;
	
	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
								
	// Read dictionary + dictionaryProps 
	const dictionary dict(particleCloud.couplingProperties());

	// Create parcel scale granular temperaturedictionary
	const dictionary propsDict(dict.subDict("ParcelScaleGranularTemperatureProps"));					
    	fileName outputRelativePath(propsDict.lookup("outputRelativePath"));
	bool verbose(false);
	if(propsDict.found("verbose")) verbose = true;			
	bool verboseParcel(false);
	if(propsDict.found("verboseParcel")) verbose = true;	
	labelList coarseGraining(propsDict.lookup("coarseGraining"));
	
        // Create output folder 
        if( !isDir(mesh.time().path()/outputRelativePath) )
        {
        	mkDir(mesh.time().path()/outputRelativePath);
        }
		
	// Create output file
	string filename;
	string prefilename = "Temp_CG_";
	std::ofstream outFile;	
	
	// OneWayVTK dict
	const dictionary VTKpropsDict(dict.subDict("oneWayVTKProps"));
	int maxNumberOfParticles=readScalar(VTKpropsDict.lookup("maxNumberOfParticles"));
	
	// Create particle cloud
	int *parcelNparts = new int[maxNumberOfParticles];
	
 	int** parcelCloud;
	
	double **parcelPositions = new double*[maxNumberOfParticles];
	for(int i = 0; i < maxNumberOfParticles; ++i) 
	{
	    parcelPositions[i] = new double[3];
	}

	double **parcelVelocities = new double*[maxNumberOfParticles];
	for(int i = 0; i < maxNumberOfParticles; ++i) 
	{
	    parcelVelocities[i] = new double[3];
	}	
		
        double **parcelKinStress = new double*[maxNumberOfParticles];
        for(int i = 0; i < maxNumberOfParticles; ++i)
        {
            parcelKinStress[i] = new double[6];
        }

	// Domain-averaged sub-parcel agitation
	scalar aveSubQparcel2(0);
					
	// Domain-averaged sub-parcel agitation
	scalar aveSuperQparcel2(0);
	
	// Mean parcel velocity
	vector 	meanParcelVel(0,0,0);

	// Total volume of domain
	scalar domainVol(0);
	forAll(mesh.C(),cellI)
	{
		domainVol +=mesh.V()[cellI];
	}
	Info << "Domain volume[m^3] = " << domainVol << endl;	
	
	// Output boolean
	// bool outputBool=false;

	// Collision algorithm

	// Dimensions, exact OR approximate  	
	int dim = 3; 

	// Number of points
	int nPts(0);
	
	// Data points
	ANNpointArray dataPts;

	// Search structure
	ANNkd_tree* kdTree;			
	
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();

		if ( runTime.timeName() != "0" )
		{
			
			int count = runTime.value() / particleCloud.dataExchangeM().DEMts();				

			Info<< " " << endl;
			particleCloud.dataExchangeM().getData("v","vector-atom",velocities,count);
			Info<< " 	Reading particle velocities" << endl;
			Info<< " " << endl;
			
			particleCloud.dataExchangeM().getData("x","vector-atom",positions,count);
			Info<< " 	Reading particle positions" << endl;
			Info<< " " << endl;
			
			particleCloud.dataExchangeM().getData("radius","scalar-atom",radii,count);
			Info<< " 	Reading particle radius" << endl;		
			Info<< " " << endl;			
		
			particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
			particleCloud.setPos(positions);
			particleCloud.setVel(velocities);

			if(verbose)
			{
				int index = 0;
				Info << "" << endl;
				Info << " index  = " << index << endl;
				Info << " rp     = " << particleCloud.radius(index) << endl;
				Info << " Vp     = " << particleCloud.velocity(index) << endl;
				Info << " Xp     = " << particleCloud.position(index) << endl;
				Info << " CellID = " << particleCloud.particleCell(index) << endl;
				Info << "" << endl;
			}
			
			
			// Number of points, collision alogrithm
			nPts =  particleCloud.numberOfParticles();	
			dataPts = annAllocPts(nPts, dim);
			
			for(int index = 0; index < particleCloud.numberOfParticles(); index++)
	        	{			
				dataPts[index][0] = particleCloud.position(index).x();		
				dataPts[index][1] = particleCloud.position(index).y();		
				dataPts[index][2] = particleCloud.position(index).z();
			}
			
			// Search structure
			kdTree = new ANNkd_tree(dataPts, nPts, dim);	
			
			for ( int parcelI = 0; parcelI < coarseGraining.size(); parcelI++ )				
			{			
				int parcelSize = coarseGraining[parcelI];
					
				// Allocate parcelCloud				
				parcelCloud = new int*[maxNumberOfParticles];
				for(int i = 0; i < maxNumberOfParticles; ++i) 
				{
	    				parcelCloud[i] = new int[parcelSize*parcelSize*parcelSize];
				}
				
				// Create parcels
				createParcels(			   mesh,
								 kdTree,
					        	  particleCloud,
							     parcelSize,
					  		    parcelCloud,
						        parcelPositions,
						       parcelVelocities,
						           parcelNparts,
							parcelKinStress,
							 aveSubQparcel2,
							  meanParcelVel,
						       	  verboseParcel		);
								
			
				// Calculate super-parcel agitation
			 	superParcelScaleTemp(     particleCloud,
						       parcelVelocities,
						          meanParcelVel,
						       aveSuperQparcel2	        );
						       
				if( verbose )
				{		       
					Info << " Coarse graining      = " << coarseGraining[parcelI];  
      					Info << " Parcel mean velocity = " << meanParcelVel;				
					Info << " aveSubQparcel2       = " << aveSubQparcel2;     
					Info << " aveSuperQparcel2     = " << aveSuperQparcel2 << endl;		       
				}
				
				filename = NumberToString(parcelSize);
				std::string fileName(mesh.time().path()/outputRelativePath/prefilename+filename);
				
				if( !outFile.is_open())
				{
					outFile.open(fileName.data(),ios_base::app);
					// if(!outputBool) outFile  << "# Time \t Mean_vel_x \t Mean_vel_y \t Mean_vel_z \t aveSubQparcel2 \t aveSuperQparcel2" << nl;				
					// outputBool=true;
				}
				
				 outFile << particleCloud.mesh().time().value()	<< "\t" 
					 << meanParcelVel[0]			<< "\t" 							
					 << meanParcelVel[1]			<< "\t" 							
					 << meanParcelVel[2]			<< "\t" 							
					 << aveSubQparcel2  			<< "\t" 		
					 << aveSuperQparcel2			<< "\t"
					 << domainVol				<< "\t"					
					 << particleCloud.numberOfParticles()	<< nl;					
								
								
				Info << " Writing into the file " << prefilename+filename << " in the folder: " << mesh.time().path()/outputRelativePath << endl;								
				
				outFile.close();
			}									
								
		}	


	}
		
	particleCloud.dataExchangeM().destroy(positions,3);
	particleCloud.dataExchangeM().destroy(velocities,3);

	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
	
	

}


// ************************************************************************* //
