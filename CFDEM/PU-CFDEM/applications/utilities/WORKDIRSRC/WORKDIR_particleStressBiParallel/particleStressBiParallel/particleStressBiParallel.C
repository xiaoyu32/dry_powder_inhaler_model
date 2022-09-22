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

#include "particleStressBiParallel.H"
#include "calcCollisionalForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	// Initiate output file
	void initOutput
	(	
		const int&          nParticleClass_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		OFstream*&          outputCollFile_,
		const bool&             bidisperse_,
		const bool&          calcCollision_		
	)
	{
	
	    	outputKinFile_ =  new OFstream(mesh_.time().path()/outputRelativePath_/"kineticStress");
		*outputKinFile_ << "#Time  " << tab << " "
		                << "Vol    " <<        " "  
	                	<< "alpp   " << tab << " "
				<< "Kin_XX " << tab << " "
				<< "Kin_YY " << tab << " "
				<< "Kin_ZZ " << tab << " " 
				<< "Kin_XY " << tab << " " 
				<< "Kin_XZ " << tab << " " 
				<< "Kin_YZ " << tab << " "; 		

		if( bidisperse_ ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        	{
				*outputKinFile_  << "Kin_XX["<<iPartClass<<"] " << tab << " "
						<< "Kin_YY["<<iPartClass<<"] " << tab << " "
						<< "Kin_ZZ["<<iPartClass<<"] " << tab << " " 
						<< "Kin_XY["<<iPartClass<<"] " << tab << " " 
						<< "Kin_XZ["<<iPartClass<<"] " << tab << " " 
						<< "Kin_YZ["<<iPartClass<<"] " << tab << " "; 
			}
		}
		*outputKinFile_ << endl;
	
		if( calcCollision_ )
		{
    			outputCollFile_ =  new OFstream(mesh_.time().path()/outputRelativePath_/"collisonalStress");
			*outputCollFile_ << "#Time  " << tab << " "
		                	 << "Vol    " << 	" "  
	                		 << "alpp   " << tab << " "			 
	                		 << "Coll_XX " << tab << " "
					 << "Coll_YY " << tab << " "
					 << "Coll_ZZ " << tab << " " 
					 << "Coll_XY " << tab << " " 
					 << "Coll_XZ " << tab << " " 
					 << "Coll_YZ " << tab << " "; 	
			if( bidisperse_ ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        		{
					*outputCollFile_  << "Coll_XX["<<iPartClass<<"] " << tab << " "
							 << "Coll_YY["<<iPartClass<<"] " << tab << " "
							 << "Coll_ZZ["<<iPartClass<<"] " << tab << " " 
							 << "Coll_XY["<<iPartClass<<"] " << tab << " " 
							 << "Coll_XZ["<<iPartClass<<"] " << tab << " " 
							 << "Coll_YZ["<<iPartClass<<"] " << tab << " "; 			
				}
			}
			*outputCollFile_ << endl;	
		}
	}
	
	// Write into output file
	void writeOutput
	(	
		const int&          nParticleClass_,
		const Time& 		   runTime_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		OFstream*&          outputCollFile_,
		const bool&             bidisperse_,
		const bool&          calcCollision_,
		SymmTensor<double>*       sigmaKin_,					
		SymmTensor<double>*	 sigmaColl_,
		const scalar&		 domainVol_,
		const scalar&	   alppInSubVolume_	)
	{
	
		// Write output
		Info << " " << endl;
		Pout << " Writing particle kinetic stresses into the file " << mesh_.time().path()/outputRelativePath_/"kineticStress" << endl;

		int iPartClass = 0;
		*outputKinFile_ << runTime_.value() 		       << tab << " " 
		                << domainVol_ 			       << tab << " "
				<< alppInSubVolume_		       << tab << " "   	
				<< sigmaKin_[iPartClass][0]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][1]/domainVol_ << tab << " " 
				<< sigmaKin_[iPartClass][2]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][3]/domainVol_ << tab << " " 
				<< sigmaKin_[iPartClass][4]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][5]/domainVol_ << tab << " " ;								

		if( bidisperse_ ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        	{
				*outputKinFile_ << sigmaKin_[iPartClass][0]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][1]/domainVol_ << tab << " " 
						<< sigmaKin_[iPartClass][2]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][3]/domainVol_ << tab << " " 
						<< sigmaKin_[iPartClass][4]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][5]/domainVol_ << tab << " " ; 
			}
		}
		*outputKinFile_ << endl;

		if( calcCollision_ )
		{
    			Pout<< " Writing particle kinetic stresses into the file " << mesh_.time().path()/outputRelativePath_/"collisionalStress" << endl;
			iPartClass = 0;
			*outputCollFile_ << runTime_.value() 		         << tab << " " 	
					 << domainVol_ 			         << tab << " "
					 << alppInSubVolume_		         << tab << " "  
					 << sigmaColl_[iPartClass][0]/domainVol_ << tab << " " 		
					 << sigmaColl_[iPartClass][1]/domainVol_ << tab << " " 
					 << sigmaColl_[iPartClass][2]/domainVol_ << tab << " " 		
					 << sigmaColl_[iPartClass][3]/domainVol_ << tab << " " 
					 << sigmaColl_[iPartClass][4]/domainVol_ << tab << " " 		
					 << sigmaColl_[iPartClass][5]/domainVol_ << tab << " " ;
			if( bidisperse_ ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        		{
					*outputCollFile_ << sigmaColl_[iPartClass][0]/domainVol_ << tab << " " 		
					  		 << sigmaColl_[iPartClass][1]/domainVol_ << tab << " " 
							 << sigmaColl_[iPartClass][2]/domainVol_ << tab << " " 		
							 << sigmaColl_[iPartClass][3]/domainVol_ << tab << " " 
							 << sigmaColl_[iPartClass][4]/domainVol_ << tab << " " 		
							 << sigmaColl_[iPartClass][5]/domainVol_ << tab << " " ;			
				}
			}
			*outputCollFile_ << endl;	
		}	
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
	double **omegas;
	double **radii;
	double **cellID;
	
	double **forces;
	double **types;

	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.dataExchangeM().allocateArray(omegas,0.,3);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
	
	particleCloud.dataExchangeM().allocateArray(forces,0.,3);
	particleCloud.dataExchangeM().allocateArray(types,0.,1);
				
	// Read dictionary + dictionaryProps 
	const dictionary dict(particleCloud.couplingProperties());
	int maxNumberOfParticles(0);
	
	// oneWayVTK dictionary
	const dictionary oneWayVTKPropsDict(dict.subDict("oneWayVTKProps"));
	maxNumberOfParticles=readScalar(oneWayVTKPropsDict.lookup("maxNumberOfParticles"));
	
	// Read particleStress Sub-dictionary
	const dictionary ParticleStressPropsDict(dict.subDict("ParticleStressProps"));
	const scalar rhop(readScalar(ParticleStressPropsDict.lookup("rhoParticle")));					
	
	// Debuging
	bool verbose(false);
	// Particle ID for debuging
	labelList exList;
	if(ParticleStressPropsDict.found("verbose")) 
	{
		verbose = true;	
		//exIndex = readInt(ParticleStressPropsDict.lookup("exIndex"));	
		exList = labelList(ParticleStressPropsDict.lookup("exList"));
	}
	
	// Collision part debugging
	bool verboseColl(false);
	if(ParticleStressPropsDict.found("verboseColl")) 
	{
		verboseColl = true;			
	}	
			
	// Number particle class
	int nParticleClass = 1;
	bool bidisperse(false);
	// Bidisperse case
	if(ParticleStressPropsDict.found("nParticleClass"))
	{
		nParticleClass = readInt(ParticleStressPropsDict.lookup("nParticleClass"));
		bidisperse = true;
		Info << " " << endl;
		Pout << "Bi-disperse case, number of particle classes = " << nParticleClass << endl;		
	}
	
	// Collision dictionary
	scalar collisionDp(1);
	bool calcCollision(false);
	if(ParticleStressPropsDict.found("calcCollision"))
	{
		calcCollision=true;
		collisionDp=readScalar(ParticleStressPropsDict.lookup("collisionDp"));
	}
	
	// Domain average
	bool domainAve(false);
	if(ParticleStressPropsDict.found("domainAve")) domainAve = true;
	

	// Parallel computation
	int nProcs(1);
	if(Pstream::parRun()) nProcs = Pstream::nProcs();

	// Create output folder	
        fileName outputRelativePath(ParticleStressPropsDict.lookup("outputRelativePath"));
	if( !isDir(mesh.time().path()/outputRelativePath) )
	{
		mkDir(mesh.time().path()/outputRelativePath );
	}

	// Initiate output file
	OFstream* outputKinFile;
	OFstream* outputCollFile;

	initOutput(	    nParticleClass,
	                              mesh,
			outputRelativePath,
			     outputKinFile,
			    outputCollFile,
			        bidisperse,
			     calcCollision	);
				
	// Define & Init collision variables
	#include "initCollisionVariables.H" 
						
	if(calcCollision)
	{
		#include "collisionParameters.H"
	}
		
	// Neighboring list parameters				
	int k(0); 

	// Dimensions, exact OR approximate  
	int dim(0); double eps(0);

	// Number of points
	int nPts(0);

	// Data points
	ANNpointArray dataPts;
	// Query points
	ANNpoint queryPt;

	ANNidxArray	nnIdx;			// 	near neighbour indices
	ANNdistArray 	dists;			//	near neighbour distances
	ANNkd_tree* 	kdTree;			//	search structure
	
	// Sub-volume of domain
	scalar domainVol(0);
	// Total volume of domain
	scalar globalDomainVol(0);	
	forAll(mesh.C(),cellI)
	{
		domainVol +=mesh.V()[cellI];
	}
	Info << " " << endl;
	Pout << "Domain volume[m^3] = " << domainVol << endl;
	globalDomainVol = domainVol;
	if(domainAve&&Pstream::parRun()) reduce(globalDomainVol, sumOp<scalar>());	
	
	// Domain Min/Max
	const pointField& pp = mesh.points();	
	
	// Min, max x-coordinates	
	scalar minX = Foam::min(pp & vector(1,0,0));
	scalar maxX = Foam::max(pp & vector(1,0,0));

	// Min, max y-coordinates		
	scalar minY = Foam::min(pp & vector(0,1,0));
	scalar maxY = Foam::max(pp & vector(0,1,0));

	// Min, max z-coordinates		
	scalar minZ = Foam::min(pp & vector(0,0,1));
	scalar maxZ = Foam::max(pp & vector(0,0,1));
	
	Pout << "x["<<minX<<":"<<maxX<<"]"<<endl;
	Pout << "y["<<minY<<":"<<maxY<<"]"<<endl;
	Pout << "z["<<minZ<<":"<<maxZ<<"]"<<endl;
	
	// Results variables
	// Collisional stress tensor
	symmTensor sigmaCollJI(0,0,0,0,0,0);
	symmTensor sigmaColl[nParticleClass+1];
	symmTensor globalSigmaColl[nParticleClass+1];
	
	// Kinetic stress tensor
	symmTensor sigmaKin[nParticleClass+1]; // nParticleClass+1 --> for mixture velocity
	symmTensor globalSigmaKin[nParticleClass+1]; 
	
	// Mean velocities in sub-volume
	vector meanVel[nParticleClass+1];	// nParticleClass+1 --> for mixture velocity
	// 
	vector globalMeanVel[nParticleClass+1];	// nParticleClass+1 --> for mixture velocity
	int npPartClass[nParticleClass+1];
	
	// Create particle list in sub-volumes
	labelList partIDInSubVolume(maxNumberOfParticles,0);

	// Create ghost particle list on boundaries
	labelList localGhostPartID(maxNumberOfParticles/nProcs+1,-1);
	labelListList ghostPartID(nProcs,localGhostPartID);
	
	// Ghost particle positions
	// Create particle cloud			
	cfdemCloud ghostParticleCloud(mesh);
	ghostParticleCloud.reAllocArrays();	
	double ** ghostPartPositions;
	ghostParticleCloud.dataExchangeM().allocateArray(ghostPartPositions,0.,3);	
		
	// Solid volume fraction in sub-volumes
	scalar alppInSubVolume;
	
	// Time loop								
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();

		if ( runTime.timeName() != "0" )
		{
			 //Read only if master
			if(Pstream::master())
			{	
				int count = runTime.value() / particleCloud.dataExchangeM().DEMts();
				dt = particleCloud.dataExchangeM().DEMts(); 				

				Info<< " " << endl;
				particleCloud.dataExchangeM().getData("v","vector-atom",velocities,count);
				Info<< " 	Reading particle velocities" << endl;
				Info<< " " << endl;

				particleCloud.dataExchangeM().getData("omega","vector-atom",omegas,count);
				Info<< " 	Reading particle angular velocities" << endl;
				Info<< " " << endl;

				particleCloud.dataExchangeM().getData("x","vector-atom",positions,count);
				Info<< " 	Reading particle positions" << endl;
				Info<< " " << endl;

				particleCloud.dataExchangeM().getData("radius","scalar-atom",radii,count);
				Info<< " 	Reading particle radius" << endl;		
				Info<< " " << endl;

				particleCloud.dataExchangeM().getData("type","scalar-atom",types,count);
				Info<< " 	Reading particle types " << endl;
				Info<< " " << endl;

				if ( verbose )
				{
					particleCloud.dataExchangeM().getData("f","vector-atom",forces,count);
					Info<< " 	Reading forces on particles" << endl;
					Info<< " " << endl;
				}
			}
			
			// Send to child CPUs
			if(Pstream::parRun())
			{
				for(int index = 0; index < maxNumberOfParticles; index++)
				{
					for(int idir = 0; idir < 3; idir++)
					{
						Pstream::scatter(velocities[index][idir]);
						Pstream::scatter(omegas[index][idir]);
						Pstream::scatter(positions[index][idir]);
					}
					Pstream::scatter(radii[index][0]);
					Pstream::scatter(types[index][0]);
				}
			}			
					
			particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
			particleCloud.setPos(positions);
			particleCloud.setVel(velocities);
			particleCloud.setOmega(omegas);
			particleCloud.setType(types);
			
			// Nunber of particles in sub-volume
			int nPartInSubVolume = 0;
			scalar totalPartVol = 0;
			int ii = 0; 
			int nGhostPart = 0;			
			
			// Search for particles touching boundaries
			for(int index = 0; index < particleCloud.numberOfParticles(); index++)
	        	{
				label cellI = particleCloud.cellIDs()[index][0];				
				if(cellI > -1)
				{
					nPartInSubVolume++;
					totalPartVol += 4./3. * constant::mathematical::pi 
							      * particleCloud.radius(index) 
							      * particleCloud.radius(index)
							      * particleCloud.radius(index); 
							      
					partIDInSubVolume[ii] = index; 
					ii++;
					
					// Find particle on boundaries
					if(       ( - minX + particleCloud.position(index).x() ) < particleCloud.radius(index) )
					{
						ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
						ghostPartPositions[nGhostPart][0] = particleCloud.position(index).x() + ( maxX - minX );
						ghostPartPositions[nGhostPart][1] = particleCloud.position(index).y() ;
						ghostPartPositions[nGhostPart][2] = particleCloud.position(index).z() ;
						
						if (verbose) Pout << " Particle " << index << " on boundary xmin " 
								  << ghostPartPositions[nGhostPart][0] << " " << particleCloud.position(index).x() << endl; 
						nGhostPart++;		  
					}
					else  if( (   maxX - particleCloud.position(index).x() ) < particleCloud.radius(index) )
					{
						ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
						ghostPartPositions[nGhostPart][0] = particleCloud.position(index).x() - ( maxX - minX );
						ghostPartPositions[nGhostPart][1] = particleCloud.position(index).y() ;
						ghostPartPositions[nGhostPart][2] = particleCloud.position(index).z() ;
						nGhostPart++;
						if (verbose) Pout << " Particle " << index << " on boundary xmax " << endl; 
					}
					else  if( ( - minY + particleCloud.position(index).y() ) < particleCloud.radius(index) )
					{
						ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
						ghostPartPositions[nGhostPart][0] = particleCloud.position(index).x() ;
						ghostPartPositions[nGhostPart][1] = particleCloud.position(index).y() + ( maxY - minY );
						ghostPartPositions[nGhostPart][2] = particleCloud.position(index).z() ;						
						nGhostPart++;
						if (verbose) Pout << " Particle " << index << " on boundary ymin " << endl; 
					}	
					else if( (   maxY - particleCloud.position(index).y() ) < particleCloud.radius(index) )
					{
						ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
						ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
						ghostPartPositions[nGhostPart][0] = particleCloud.position(index).x() ;
						ghostPartPositions[nGhostPart][1] = particleCloud.position(index).y() - ( maxY - minY );
						ghostPartPositions[nGhostPart][2] = particleCloud.position(index).z() ;
						nGhostPart++;
						if (verbose) Pout << " Particle " << index << " on boundary ymax " << endl;
					}	 
					else if( ( - minZ + particleCloud.position(index).z() ) < particleCloud.radius(index) )
					{
						ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
						ghostPartPositions[nGhostPart][0] = particleCloud.position(index).x() ;
						ghostPartPositions[nGhostPart][1] = particleCloud.position(index).y() ;
						ghostPartPositions[nGhostPart][2] = particleCloud.position(index).z() + ( maxZ - minZ );					
						nGhostPart++;
						if (verbose) Pout << " Particle " << index << " on boundary zmin " << endl;
					}	 
					else if( (   maxZ - particleCloud.position(index).z() ) < particleCloud.radius(index) )
					{
						ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
						ghostPartPositions[nGhostPart][0] = particleCloud.position(index).x() ;
						ghostPartPositions[nGhostPart][1] = particleCloud.position(index).y() ;
						ghostPartPositions[nGhostPart][2] = particleCloud.position(index).z() - ( maxZ - minZ );
						nGhostPart++;
						if (verbose) Pout << " Particle " << index << " on boundary zmax " << endl; 						
					}		            							      
				}
			}	
			
			// Resize partID list in sub-volume
			partIDInSubVolume.resize(nPartInSubVolume);
			alppInSubVolume = totalPartVol / domainVol;
			if(Pstream::parRun()) Pout << " Number of particles in sub-volume = " << nPartInSubVolume <<endl;
			if(Pstream::parRun()) Pout << " Solid volume fraction in sub-volume = " << alppInSubVolume <<endl;
			
			// Parallel computation
			// Resize number of ghost particle array
			ghostPartID[Pstream::myProcNo()].resize(nGhostPart);
									
			// Distrubute to all processors
			Pstream::gatherList(ghostPartID);
			Pstream::scatterList(ghostPartID);
			
			// Create global ghost particle IDs
			labelList globalGhostPartID(nGhostPart,-1);
			
			globalGhostPartID = ListListOps::combine<labelList>
			(
				ghostPartID,
				accessOp<labelList>()
			);
			
			// Set number of parcels
			ghostParticleCloud.setNumberOfParticles(globalGhostPartID.size());	 

			/*
			if(Pstream::parRun())
			{			
				for(int index = 0; index < globalGhostPartID.size(); index++)
				{
					for(int idir=0; idir<3; idir++)
					{
						Pstream::mapCombineScatter(ghostPartPositions[index][idir]);
					}
				}	
			}
			*/
					
			
			// Set ghost particle positions
			//ghostParticleCloud.setPos(ghostPartPositions);
			//ghostParticleCloud.setPos(fcoll);
				
			//
			if(verbose)
			{
				Pout << " Ghost ID " << tab 
				     << "Position[x]" << tab << "Real pos.[x]" << tab
				     << "Position[y]" << tab << "Real pos.[y]" << tab
				     << "Position[z]" << tab << "Real pos.[z]" << endl;;
				for(int index = 0; index < globalGhostPartID.size(); index++)
				{
					Pout << index << tab;
					Pout << globalGhostPartID[index] << tab;
					for(int idir=0; idir<3; idir++)
					{
						//Pout << ghostParticleCloud.position(index)[idir] << tab;
						Pout << ghostPartPositions[index][idir] << tab;
						Pout << particleCloud.position(globalGhostPartID[index])[idir] <<tab;
					}
					Pout << endl;
				}	
			}
			
			if(verbose)
			{
				for(label ii=0; ii < exList.size(); ii++)
				{
					Info << "" << endl;
					int exIndex = exList[ii]; 
					label cellI = particleCloud.cellIDs()[exIndex][0];				
					if(cellI > -1)
					{					
						Pout << " index  = " << exIndex << endl;
						Pout << " rp     = " << particleCloud.radius(exIndex) << endl;
						Pout << " Type   = " << particleCloud.type(exIndex) << endl;
						Pout << " Vp     = " << particleCloud.velocity(exIndex) << endl;
						Pout << " Omegap = " << particleCloud.omega(exIndex) << endl;				
						Pout << " Xp     = " << particleCloud.position(exIndex) << endl;
						Pout << " Fp     = " << forces[exIndex][0] << " " << forces[exIndex][1] << " " << forces[exIndex][2] << endl;
					}	
				}
			}				
			
			if(calcCollision)
			{
				// Neighboring list parameters
				k = 20; 
				//if(particleCloud.numberOfParticles()<k) k = particleCloud.numberOfParticles();
				if(nPartInSubVolume < k) k = nPartInSubVolume; 

				// Dimensions, exact OR approximate  
				dim=3; eps = 0;

				// Number of points
				//nPts =  particleCloud.numberOfParticles();		       		
				nPts = nPartInSubVolume;  

				// Allocate 
				queryPt = annAllocPt(dim);
				dataPts = annAllocPts(nPts, dim);
				nnIdx = new ANNidx[k];
				dists = new ANNdist[k];

				// Particle collisional stresses
				sigmaCollJI = symmTensor(0,0,0,0,0,0);				
				for(int iPartClass = 0; iPartClass <= nParticleClass; iPartClass++)
	        		{
					sigmaColl[iPartClass] = symmTensor(0,0,0,0,0,0);
					if(domainAve) globalSigmaColl[iPartClass] = symmTensor(0,0,0,0,0,0);
				}	
			}			
			
			
						
			// Initiate
			int iPartClass = 0;
			//npPartClass[iPartClass] = particleCloud.numberOfParticles();
			npPartClass[iPartClass] = nPartInSubVolume;  
			
			sigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0);
			globalSigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0); 
			
			if( bidisperse ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass; iPartClass++)
	        		{
					// Initiate particle number of each class
					npPartClass[iPartClass] = 0;
					
					// Initiate velocities
					for(int idir=0; idir<3; idir++) 
					{
						meanVel[iPartClass][idir] = 0 ;
						if(domainAve) globalMeanVel[iPartClass][idir] = 0 ;
					}

					// Initiate particle kinetic stresses
					sigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0);
					// Initiate global particle kinetic stresses
					globalSigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0);
				}
			}		
				
			// Create particle list
			//for(int index = 0; index <  particleCloud.numberOfParticles(); index++)
			for(int ii = 0; ii <  partIDInSubVolume.size(); ii++)
			{							
				   label index = partIDInSubVolume[ii];	
				//  Cell ID
				// label cellI = particleCloud.cellIDs()[index][0];
				// if(cellI > -1)
				// {
					if(calcCollision)	
					{
						//dataPts[index][0] = particleCloud.position(index).x();		
						//dataPts[index][1] = particleCloud.position(index).y();		
						//dataPts[index][2] = particleCloud.position(index).z();
						
						dataPts[ii][0] = particleCloud.position(index).x();		
						dataPts[ii][1] = particleCloud.position(index).y();		
						dataPts[ii][2] = particleCloud.position(index).z();						

						for (int dir=0;dir<3;dir++)
						{	
							fcoll[index][dir] = 0;
					 		ftan[index][dir] = 0;		
 					 		fcap[index][dir] = 0;
 							fvisc[index][dir] = 0;
						}
					}

					// Total velocity of particles
					int iPartClass = 0;
					for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += particleCloud.velocity(index)[idir];

					if( bidisperse ) 
					{
						iPartClass = particleCloud.type(index);

						for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += particleCloud.velocity(index)[idir];
						npPartClass[iPartClass]++;	
					}	
				//}
			}	
	
			// Global domain average
			if (domainAve) 
			{
				for(int iPartClass = 0; iPartClass <= nParticleClass ; iPartClass++)
	        		{
					for(int idir=0; idir<3; idir++) globalMeanVel[iPartClass][idir] = meanVel[iPartClass][idir];

					// Parallel computation
					reduce(globalMeanVel[iPartClass], sumOp<vector>());
					// Normalize
					for(int idir=0; idir<3; idir++) globalMeanVel[iPartClass][idir] /= particleCloud.numberOfParticles();
				}
			}			
			
			// Normalize
			for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir]/= npPartClass[iPartClass];
			if( bidisperse ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	        		{
					for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir]/= npPartClass[iPartClass];
				}
			}				
			if(verbose)
			{
				iPartClass = 0;									
				Info << " " << endl;
				Pout << " Particle class = " << iPartClass << endl;
				Pout << " <u_p,x> = " << meanVel[iPartClass][0] << endl;
				Pout << " <u_p,y> = " << meanVel[iPartClass][1] << endl;
				Pout << " <u_p,z> = " << meanVel[iPartClass][2] << endl;

				if (domainAve && Pstream::master())
				{				
					Pout << " Domain <u_p,x> = " << globalMeanVel[iPartClass][0] << endl;
					Pout << " Domain <u_p,y> = " << globalMeanVel[iPartClass][1] << endl;
					Pout << " Domain <u_p,z> = " << globalMeanVel[iPartClass][2] << endl;
				}
				
				if( bidisperse ) 
				{
					for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	        			{
						Info << " " << endl;
						Pout << " Particle class = " << iPartClass << endl;
						Pout << " <u_p,x> = " << meanVel[iPartClass][0] << endl;
						Pout << " <u_p,y> = " << meanVel[iPartClass][1] << endl;
						Pout << " <u_p,z> = " << meanVel[iPartClass][2] << endl;
						
						if (domainAve && Pstream::master())
						{				
							Pout << " Domain <u_p,x> = " << globalMeanVel[iPartClass][0] << endl;
							Pout << " Domain <u_p,y> = " << globalMeanVel[iPartClass][1] << endl;
							Pout << " Domain <u_p,z> = " << globalMeanVel[iPartClass][2] << endl;
						}
												
					}
				}
								
			}

			if(calcCollision)
			{
				// Create Tree structure
				kdTree = new ANNkd_tree(dataPts, nPts, dim);
			}
			
			vector pos_i(0,0,0);
			vector pos_j(0,0,0);
			scalar dist(0);
			
			// Particle mass
			scalar mass_p(0);
			int typeJ(0);
			int typeI(0);
			label cellI;
					
			//for(int index_j = 0; index_j <  particleCloud.numberOfParticles(); index_j++)		
			for(int ii = 0; ii <  partIDInSubVolume.size(); ii++)
			{
				// Cell ID
				//cellI = particleCloud.cellIDs()[index_j][0];
				label index_j = partIDInSubVolume[ii]; 
				
				// Particle type
				typeJ = particleCloud.type(index_j) - 1; // Just nclass starts from "0"
				
				//if(cellI > -1)
				//{
					if(calcCollision)
					{
						scalar sqRad = collisionDp*particleCloud.radius(index_j);

						queryPt[0] = particleCloud.position(index_j).x();
						queryPt[1] = particleCloud.position(index_j).y();
						queryPt[2] = particleCloud.position(index_j).z();

						kdTree->annkFRSearch(
			                				queryPt,			// query point					
									sqRad,				// squared radius
									k,				// number of the near neighbours to return
									nnIdx,				// nearest neighbor array
									dists,				// dist to near neighbours
									eps			);					

						//for (int index_i = 1; index_i < k; index_i++) // k=0 particle itself
						for (int index_i = 0; index_i < k; index_i++) // k=0 particle itself ???
						{

							if ( nnIdx[index_i] != index_j )
							{
						         typeI = particleCloud.type(nnIdx[index_i])-1; // Just nclass stats from "0"

									// Calculate collision forces
									calcForce(			  
											particleCloud,
										      collisionModelI, 
											      index_j,
										       nnIdx[index_i],     
												fcoll, 
												 ftan,
												  k_n, 
												  k_t, 
											      gamma_n, 
											      gamma_t, 
                                                        				youngsModulus,
                                                        				poissonsRatio,
                                        				       coefficientRestitution,
                                                				  coefficientFriction,				
											      delta_t, 
												 mu_f, 
												   dt, 
					        				   tangential_history, 
							  					  liq, 
										      liquid_transfer, 
						        				    liquidVol, 
						        				 surf_tension, 
						        				   fluid_visc, 
							        				 fcap, 
												fvisc,
						        				  first_touch,
											     cohesion,
										   minimumDistanceVdW,
						        				cohEnergyDens,
												 fcoh,
												 rhop,
												  e_n,
												  e_t,
											  sigmaCollJI,
											   bidisperse,
								        			typeJ,
												typeI,
											  verboseColl   		);


									iPartClass = 0;
                                                			sigmaColl[iPartClass][0] += sigmaCollJI.xx();
                                                			sigmaColl[iPartClass][1] += sigmaCollJI.yy();
                                                			sigmaColl[iPartClass][2] += sigmaCollJI.zz();
                                                			sigmaColl[iPartClass][3] += sigmaCollJI.xy();
                                                			sigmaColl[iPartClass][4] += sigmaCollJI.xz();
                                                			sigmaColl[iPartClass][5] += sigmaCollJI.yz();


									if( bidisperse )
									{							  						
                                                        			if( typeJ == typeI )
                                                        			{
                                                 
											iPartClass = typeJ + 1 ;
											sigmaColl[iPartClass][0] += sigmaCollJI.xx();
                                                                			sigmaColl[iPartClass][1] += sigmaCollJI.yy();
                                                                			sigmaColl[iPartClass][2] += sigmaCollJI.zz();
                                                                			sigmaColl[iPartClass][3] += sigmaCollJI.xy();
                                                                			sigmaColl[iPartClass][4] += sigmaCollJI.xz();
                                                                			sigmaColl[iPartClass][5] += sigmaCollJI.yz();
                                                        			}

									}
									
							}
						
						}	
							
						// Loop over ghost particles 					
						for (int index_i = 0; index_i < globalGhostPartID.size(); index_i++) //
						{

							if ( globalGhostPartID[index_i] != index_j )
							{
						         typeI = particleCloud.type(globalGhostPartID[index_i])-1; // Just nclass stats from "0"
									
									// Calculate collision forces
									calcForce(			  
											particleCloud,
										      collisionModelI, 
											      index_j,
								           globalGhostPartID[index_i],     
												fcoll, 
												 ftan,
												  k_n, 
												  k_t, 
											      gamma_n, 
											      gamma_t, 
                                                        				youngsModulus,
                                                        				poissonsRatio,
                                        				       coefficientRestitution,
                                                				  coefficientFriction,				
											      delta_t, 
												 mu_f, 
												   dt, 
					        				   tangential_history, 
							  					  liq, 
										      liquid_transfer, 
						        				    liquidVol, 
						        				 surf_tension, 
						        				   fluid_visc, 
							        				 fcap, 
												fvisc,
						        				  first_touch,
											     cohesion,
										   minimumDistanceVdW,
						        				cohEnergyDens,
												 fcoh,
												 rhop,
												  e_n,
												  e_t,
											  sigmaCollJI,
											   bidisperse,
								        			typeJ,
												typeI,
											  verboseColl   		);


									iPartClass = 0;
                                                			sigmaColl[iPartClass][0] += sigmaCollJI.xx();
                                                			sigmaColl[iPartClass][1] += sigmaCollJI.yy();
                                                			sigmaColl[iPartClass][2] += sigmaCollJI.zz();
                                                			sigmaColl[iPartClass][3] += sigmaCollJI.xy();
                                                			sigmaColl[iPartClass][4] += sigmaCollJI.xz();
                                                			sigmaColl[iPartClass][5] += sigmaCollJI.yz();


									if( bidisperse )
									{							  						
                                                        			if( typeJ == typeI )
                                                        			{
                                                 
											iPartClass = typeJ + 1 ;
											sigmaColl[iPartClass][0] += sigmaCollJI.xx();
                                                                			sigmaColl[iPartClass][1] += sigmaCollJI.yy();
                                                                			sigmaColl[iPartClass][2] += sigmaCollJI.zz();
                                                                			sigmaColl[iPartClass][3] += sigmaCollJI.xy();
                                                                			sigmaColl[iPartClass][4] += sigmaCollJI.xz();
                                                                			sigmaColl[iPartClass][5] += sigmaCollJI.yz();
                                                        			}

									}
									
							}
					
							
								

						}

					}	

					// Particle mass
					mass_p = 4./3.*rhop*constant::mathematical::pi*particleCloud.radius(index_j)*particleCloud.radius(index_j)*particleCloud.radius(index_j);	

					// Particle kinetic stress
					int iPartClass = 0;
					sigmaKin[iPartClass][0] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                        	  * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ); 

					sigmaKin[iPartClass][1] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
									  * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ); 

					sigmaKin[iPartClass][2] += mass_p * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ) 
									  * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ); 

					sigmaKin[iPartClass][3] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                        	  * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] );

					sigmaKin[iPartClass][4] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                        	  * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );

					sigmaKin[iPartClass][5] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
									  * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );
									  
					if(domainAve)
					{
						globalSigmaKin[iPartClass][0] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
					                        		        * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ); 

						globalSigmaKin[iPartClass][1] += mass_p * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ) 
										        * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ); 

						globalSigmaKin[iPartClass][2] += mass_p * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] ) 
										        * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] ); 

						globalSigmaKin[iPartClass][3] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
					                        		        * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] );

						globalSigmaKin[iPartClass][4] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
					                        		        * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] );

						globalSigmaKin[iPartClass][5] += mass_p * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ) 
										        * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] );						
					
					}				  
											   
					if( bidisperse ) 
					{				
						iPartClass = typeJ + 1;

						sigmaKin[iPartClass][0] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                                	  * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ); 

						sigmaKin[iPartClass][1] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
										  * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ); 

						sigmaKin[iPartClass][2] += mass_p * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ) 
										  * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ); 

						sigmaKin[iPartClass][3] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                                	  * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] );

						sigmaKin[iPartClass][4] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                                	  * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );

						sigmaKin[iPartClass][5] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
										  * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );	
										  
						if(domainAve)
						{
							globalSigmaKin[iPartClass][0] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
					                        		        	* ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ); 

							globalSigmaKin[iPartClass][1] += mass_p * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ) 
										        	* ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ); 

							globalSigmaKin[iPartClass][2] += mass_p * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] ) 
										        	* ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] ); 

							globalSigmaKin[iPartClass][3] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
					                        		        	* ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] );

							globalSigmaKin[iPartClass][4] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
					                        		        	* ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] );

							globalSigmaKin[iPartClass][5] += mass_p * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ) 
										        	* ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] );						

						}
										   
					}	
				//}
			}
			
			
			// Write output
			writeOutput(	    nParticleClass,
						   runTime,
			                              mesh,
					outputRelativePath,
			     		     outputKinFile,
			    		    outputCollFile,
			        		bidisperse,
			     		     calcCollision,
					          sigmaKin,
					         sigmaColl,
						 domainVol,
					   alppInSubVolume	);		   
						 
			if( domainAve )
			{
				for(int iPartClass = 0; iPartClass <= nParticleClass ; iPartClass++)
	        		{
					for(label ii=0; ii < 6; ii++)
					{
						globalSigmaColl[iPartClass][ii] = sigmaColl[iPartClass][ii];
						// Parallel computation
						reduce(globalSigmaKin[iPartClass][ii], sumOp<scalar>());
						reduce(globalSigmaColl[iPartClass][ii], sumOp<scalar>());					
					}
				}
			}

			if(verbose)
			{
				if(calcCollision)
				{
					for(label ii=0; ii < exList.size(); ii++)
					{
						Info << "" << endl;
						int exIndex = exList[ii]; 
						Info << "" << endl;
						if ( particleCloud.cellIDs()[exIndex][0] > -1 )
				 			Info << " Fcoll = " << fcoll[exIndex][0] << " " << fcoll[exIndex][1] << " " << fcoll[exIndex][2] << endl;
					}
				}
				if(cohesion)
				{
					for(label ii=0; ii < exList.size(); ii++)
					{
						Info << "" << endl;
						int exIndex = exList[ii]; 
						if ( particleCloud.cellIDs()[exIndex][0] > -1 )
							Info << " Fcoh = " <<  fcoh[exIndex][0] << " " <<  fcoh[exIndex][1] << " " <<  fcoh[exIndex][2] << endl;
						Info << "" << endl;
					}	
				}	
			
				iPartClass = 0;
				Info << " " << endl;
				Pout << " Particle class = " << iPartClass << " kinetic stresses" << endl;
				Pout << " sigmaKin_xx= " << sigmaKin[iPartClass][0]/domainVol << endl;
				Pout << " sigmaKin_yy= " << sigmaKin[iPartClass][1]/domainVol << endl;
				Pout << " sigmaKin_zz= " << sigmaKin[iPartClass][2]/domainVol << endl;			
				Pout << " sigmaKin_xy= " << sigmaKin[iPartClass][3]/domainVol << endl;
				Pout << " sigmaKin_xz= " << sigmaKin[iPartClass][4]/domainVol << endl;
				Pout << " sigmaKin_yz= " << sigmaKin[iPartClass][5]/domainVol << endl;

				if( domainAve && Pstream::master() )
				{						
					Info << " Domain sigmaKin_xx= " << globalSigmaKin[iPartClass][0]/globalDomainVol << endl;
					Info << " Domain sigmaKin_yy= " << globalSigmaKin[iPartClass][1]/globalDomainVol << endl;
					Info << " Domain sigmaKin_zz= " << globalSigmaKin[iPartClass][2]/globalDomainVol << endl;			
					Info << " Domain sigmaKin_xy= " << globalSigmaKin[iPartClass][3]/globalDomainVol << endl;
					Info << " Domain sigmaKin_xz= " << globalSigmaKin[iPartClass][4]/globalDomainVol << endl;
					Info << " Domain sigmaKin_yz= " << globalSigmaKin[iPartClass][5]/globalDomainVol << endl;													
				}
				
				if(calcCollision)
				{
					Info << " " << endl;
					Pout << " Particle class = " << iPartClass << " collisional stresses" << endl;
					Pout << " sigmaColl_xx= " << sigmaColl[iPartClass][0]/domainVol << endl;
					Pout << " sigmaColl_yy= " << sigmaColl[iPartClass][1]/domainVol << endl;
					Pout << " sigmaColl_zz= " << sigmaColl[iPartClass][2]/domainVol << endl;						
					Pout << " sigmaColl_xy= " << sigmaColl[iPartClass][3]/domainVol << endl;
					Pout << " sigmaColl_xz= " << sigmaColl[iPartClass][4]/domainVol << endl;			
					Pout << " sigmaColl_yz= " << sigmaColl[iPartClass][5]/domainVol << endl;
					
					if( domainAve && Pstream::master() )
					{						
						Info << " Domain sigmaColl_xx= " << globalSigmaColl[iPartClass][0]/globalDomainVol << endl;
						Info << " Domain sigmaColl_yy= " << globalSigmaColl[iPartClass][1]/globalDomainVol << endl;
						Info << " Domain sigmaColl_zz= " << globalSigmaColl[iPartClass][2]/globalDomainVol << endl;						
						Info << " Domain sigmaColl_xy= " << globalSigmaColl[iPartClass][3]/globalDomainVol << endl;
						Info << " Domain sigmaColl_xz= " << globalSigmaColl[iPartClass][4]/globalDomainVol << endl;			
						Info << " Domain sigmaColl_yz= " << globalSigmaColl[iPartClass][5]/globalDomainVol << endl;														
					}

				}
				

				
				if( bidisperse ) 
				{
					for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	        			{
						Info << " " << endl;
						Pout << " Particle class = " << iPartClass << " kinetic stresses" << endl;
						Pout << " sigmaKin_xx= " << sigmaKin[iPartClass][0]/domainVol << endl;
						Pout << " sigmaKin_yy= " << sigmaKin[iPartClass][1]/domainVol << endl;
						Pout << " sigmaKin_zz= " << sigmaKin[iPartClass][2]/domainVol << endl;			
						Pout << " sigmaKin_xy= " << sigmaKin[iPartClass][3]/domainVol << endl;
						Pout << " sigmaKin_xz= " << sigmaKin[iPartClass][4]/domainVol << endl;
						Pout << " sigmaKin_yz= " << sigmaKin[iPartClass][5]/domainVol << endl;
						
						if( domainAve && Pstream::master() )
						{						
							Info << " Domain sigmaKin_xx= " << globalSigmaKin[iPartClass][0]/globalDomainVol << endl;
							Info << " Domain sigmaKin_yy= " << globalSigmaKin[iPartClass][1]/globalDomainVol << endl;
							Info << " Domain sigmaKin_zz= " << globalSigmaKin[iPartClass][2]/globalDomainVol << endl;			
							Info << " Domain sigmaKin_xy= " << globalSigmaKin[iPartClass][3]/globalDomainVol << endl;
							Info << " Domain sigmaKin_xz= " << globalSigmaKin[iPartClass][4]/globalDomainVol << endl;
							Info << " Domain sigmaKin_yz= " << globalSigmaKin[iPartClass][5]/globalDomainVol << endl;												
						}						

						if(calcCollision)
						{
							Info << " " << endl;
							Pout << " Particle class = " << iPartClass << " collisional stresses" << endl;
							Pout << " sigmaColl_xx= " << sigmaColl[iPartClass][0]/domainVol << endl;
							Pout << " sigmaColl_yy= " << sigmaColl[iPartClass][1]/domainVol << endl;
							Pout << " sigmaColl_zz= " << sigmaColl[iPartClass][2]/domainVol << endl;						
							Pout << " sigmaColl_xy= " << sigmaColl[iPartClass][3]/domainVol << endl;
							Pout << " sigmaColl_xz= " << sigmaColl[iPartClass][4]/domainVol << endl;			
							Pout << " sigmaColl_yz= " << sigmaColl[iPartClass][5]/domainVol << endl;
							
							if( domainAve && Pstream::master() )
							{						
								Info << " Domain sigmaColl_xx= " << globalSigmaColl[iPartClass][0]/globalDomainVol << endl;
								Info << " Domain sigmaColl_yy= " << globalSigmaColl[iPartClass][1]/globalDomainVol << endl;
								Info << " Domain sigmaColl_zz= " << globalSigmaColl[iPartClass][2]/globalDomainVol << endl;						
								Info << " Domain sigmaColl_xy= " << globalSigmaColl[iPartClass][3]/globalDomainVol << endl;
								Info << " Domain sigmaColl_xz= " << globalSigmaColl[iPartClass][4]/globalDomainVol << endl;			
								Info << " Domain sigmaColl_yz= " << globalSigmaColl[iPartClass][5]/globalDomainVol << endl;														
							}

						}						
					}
				}


			
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
