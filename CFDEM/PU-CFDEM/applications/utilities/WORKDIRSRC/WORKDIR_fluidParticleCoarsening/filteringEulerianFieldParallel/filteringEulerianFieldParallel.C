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

#include "filteringEulerianFieldParallel.H"
#include "domainFilter.H"
#include "readVariables.H"
#include "dragForce.H"
#include "fluidCoarsening.H"
#include "createParcel.H"
#include "coarsening.H"
#include "applyBin.H"
#include "output.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
		
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
	
	Info << " Preparing the case..." << endl;

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
	
	
	// Create parcel cloud	( Max:5 )
	cfdemCloud parcelCloud(mesh);
	parcelCloud.reAllocArrays();
	
	/*
	cfdemCloud parcelCloud_2(mesh);
	parcelCloud_2.reAllocArrays();
	
	cfdemCloud parcelCloud_3(mesh);
	parcelCloud_3.reAllocArrays();
	
	cfdemCloud parcelCloud_4(mesh);
	parcelCloud_4.reAllocArrays();
	
	cfdemCloud parcelCloud_5(mesh);
	parcelCloud_5.reAllocArrays();					
	*/
		
	// Parcel velocities are defined in parcel coarsening subroutine	
	
	// Post-processing filtered drag dictionary 
	// Read dictionary + dictionaryProps 
	const dictionary dict(particleCloud.couplingProperties());
	// Data exchange oneWayVTK
	/*const dictionary VTKDict(dict.subDict("oneWayVTKProps"));					
	const scalar DEMts(readScalar(VTKDict.lookup("DEMts")));
    	fileName relativePath(VTKDict.lookup("relativePath"));
	fileName couplingFilename(VTKDict.lookup("couplingFilename"));
	const int maxNumberOfParticles(readScalar(VTKDict.lookup("maxNumberOfParticles")));	*/
				
	// Wen&Yu Drag
	const dictionary propsDict(dict.subDict("PostProcessingFilteredDrag"));					
	bool verbose(false);
	//const int DEM_dump_Interval(readScalar(propsDict.lookup("DEM_dump_Interval")));	
	const int nBin(readScalar(propsDict.lookup("nBin")));
	const scalar maxalpp(readScalar(propsDict.lookup("Max_alpp")));
	const scalar minalpp(readScalar(propsDict.lookup("Min_alpp")));
	// Binning with VelSlip, taup |Vr|
	//const scalar minVelSlip(-1);
	//const scalar maxVelSlip( 1);
	//const scalar minTaupVelSlip(0);
	//const scalar maxTaupVelSlip(0.5);
	//
	const scalar minVelSlip(readScalar(propsDict.lookup("Min_VelSlip")));
	const scalar maxVelSlip(readScalar(propsDict.lookup("Max_VelSlip")));
	const scalar minTaupVelSlip(readScalar(propsDict.lookup("Min_TaupVelSlip")));
	const scalar maxTaupVelSlip(readScalar(propsDict.lookup("Max_TaupVelSlip")));	
	//
	const scalar rhop(readScalar(propsDict.lookup("rhopart")));
    	fileName outputRelativePath(propsDict.lookup("outputRelativePath"));
	fileName stencilsOutputRelativePath(propsDict.lookup("stencilsOutputRelativePath"));
	bool particleCoarsening(false);
	bool verboseParticleCoarsening(false);
	bool weighting(false);
	//int npartParticleCoarsening(1);
	labelList npartParticleCoarsening(1);
	if(propsDict.found("verbose")) verbose = true;	
	if(propsDict.found("ParticleCoarsening")) 
	{
		particleCoarsening = true;
		if(propsDict.found("verboseParticleCoarsening")) verboseParticleCoarsening = true; 
		npartParticleCoarsening = propsDict.lookup("npartParticleCoarsening");
		labelList test(propsDict.lookup("npartParticleCoarsening"));
		if ( test.size() > 5 )  FatalError << " Parcel list size cannot be greater than 5 " << abort(FatalError);
	}
	
	if(propsDict.found("weighting")) weighting = true;
	
	// Filtering Euler-Euler Approach
	bool EulerEulerFiltering(false);
	if(propsDict.found("EulerEulerFiltering")) 
	{
		Info << "\nEuler-Euler filtering is active " << endl;
		EulerEulerFiltering = true;
	}
	
	// Euler-Lagrange filtering with the second binning Vslip in Eulerian manner
	bool EulerianVslipBin(false);
	if(propsDict.found("EulerianVslipBin")) 
	{
		Info << "\nThe second binning Vslip in Eulerian manner " << endl;
		EulerianVslipBin = true;
	}
		
	// Filter parameters
	#include "FilterVariables.H"	

	if(!Pstream::parRun())
	{	
		Foam::constructfilter 
		(
			args,
			runTime,
			mesh,
			minFilterWidth,
			maxFilterWidth,
			FilterIncrement,
			StencilListFilter,
			stencilsOutputRelativePath
		);
	}						
	
	// ######################### NEW VARIABLES (Parallel compuatation) #######################

	// oneWayVTK dictionary
	const dictionary oneWayVTKPropsDict(dict.subDict("oneWayVTKProps"));
	int maxNumberOfParticles=readScalar(oneWayVTKPropsDict.lookup("maxNumberOfParticles"));

	// Domain average
	bool domainAve(false);
	if(propsDict.found("domainAve")) domainAve = true;
		
	// Create particle list in sub-volumes
	labelList partIDInSubVolume(maxNumberOfParticles,-1);

	// Create parcel list in sub-volumes
	labelList parcelIDInSubVolume(maxNumberOfParticles,-1);
	
	// Solid volume fraction in sub-volumes
	scalar alppInSubVolume(0);

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
	
	// ######################### NEW VARIABLES (Parallel compuatatio) #######################
					
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();

		if ( runTime.timeName() != "0" )
		{
			readEulerianVariables
			(
				args, 
				runTime, 
				mesh,
				voidfraction,
				U,
				rho,
				p
			);
			
			 //Read only if master
			if(Pstream::master())
			{
				int count = runTime.value() / particleCloud.dataExchangeM().DEMts();				
				// timeI > 0 (particles are not initialised in the folder "0")

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
			}	
			
			// Send to child CPUs
			if(Pstream::parRun())
			{
				for(int index = 0; index < maxNumberOfParticles; index++)
				//for(int index = 0; index < particleCloud.numberOfParticles(); index++)
				{
					for(int idir = 0; idir < 3; idir++)
					{
						Pstream::scatter(velocities[index][idir]);
						Pstream::scatter(positions[index][idir]);
					}
					Pstream::scatter(radii[index][0]);
				}
			}
		
			particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
			particleCloud.setPos(positions);
			particleCloud.setVel(velocities);

			// Nunber of particles in sub-volume
			int nPartInSubVolume = 0;
			scalar totalPartVol = 0;
			int ii = 0; 
						
			// Particle IDs in sub-volume
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
					// Pout << ii <<  " cellI " << cellI  << " partIDInSubVolume[ii] "	<< partIDInSubVolume[ii] << endl;         					
					ii++;				
				}
			}	
					
			// Resize partID list in sub-volume
			partIDInSubVolume.resize(nPartInSubVolume);				
			alppInSubVolume = totalPartVol / domainVol;
			if(Pstream::parRun()) Pout << " Number of particles in sub-volume = " << nPartInSubVolume <<endl;
			if(Pstream::parRun()) Pout << " Solid volume fraction in sub-volume = " << alppInSubVolume <<endl;

			// Create particle list
			labelList particleL(1,(0));
			labelListList particleList(partIDInSubVolume.size(),particleL);
			for(int ii = 0; ii < partIDInSubVolume.size(); ii++)
			{
				particleList[ii][0] = partIDInSubVolume[ii];
			}		

			if(verbose)
			{
				int index = 0;
				if( particleCloud.cellIDs()[index][0] > -1 )
				{	
					Info << "" << endl;
					Pout << " index  = " << index << endl;
					Pout << " rp     = " << particleCloud.radius(index) << endl;
					Pout << " Vp     = " << particleCloud.velocity(index) << endl;
					Pout << " Xp     = " << particleCloud.position(index) << endl;
					Pout << " CellID = " << particleCloud.particleCell(index) << endl;
					Pout << " Ug     = " << U[particleCloud.particleCell(index)] << endl;
					Info << "" << endl;
				}	
			}

			/* Check cellID numbering!!! local numberinghas to be called as particleCloud.cellIDs()[ partIDInSubVolume[ii]][0]
			for(int ii =0; ii < partIDInSubVolume.size(); ii++)
			{
				Pout << " Particle = " << ii << " Global particle ID " << partIDInSubVolume[ii]
				     << " CellID " << particleCloud.cellIDs()[ii][0]  << " Global CellID " << particleCloud.cellIDs()[ partIDInSubVolume[ii]][0] << endl;
			}	
			*/			
			
			// 	Call DragForce function for the drag force	
			CalculateDragForce
			(
				particleCloud,
				voidfraction,
				U,
				rho,
				verbose,
				DragForce,
				partIDInSubVolume
			);				
						
			if( EulerEulerFiltering || EulerianVslipBin ) 
			{
				Info << " CPCCellToCellStencil not working in parallel computation... " << endl;
									
				Foam::EulerianParticleVelocityForce
				(
					particleCloud,			
					mesh,
					Up,
					DragForce,
					MappedDragForce,
					partIDInSubVolume			
				);								
								
				Up.write();								
					
				Foam::CalculateEulerianDragForce
				(
					particleCloud,
					mesh,
					voidfraction,
					U,
					Up,
					rho,
					EulerianDragForce	
				);

			}

			// Parcel particle list
			labelListListList parcelParticleListList;
			
			if(particleCoarsening)
			{
				//forAll(npartParticleCoarsening,npartI)
				for ( int npartI = 0; npartI < npartParticleCoarsening.size(); npartI++ )				
				{
					Info << " " << endl;
					Info << " Particle coarsening: "<< endl;
					
					// Number of parcel
					//int nParcel =  particleCloud.numberOfParticles()/npartParticleCoarsening[npartI];					
					int nParcel =  partIDInSubVolume.size()/npartParticleCoarsening[npartI];	
					labelList particleID(npartParticleCoarsening[npartI],(0));
					labelListList parcelParticleList(nParcel,particleID);
					//parcelParticleListList.setSize(5,parcelParticleList);		
					parcelParticleListList.setSize(npartParticleCoarsening.size(),parcelParticleList);
					
					Info << " " << endl;				
					Pout << "Creating parcels..." << endl;
					createParcels(mesh,particleCloud,parcelCloud,npartParticleCoarsening[npartI],
					  			  verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting,
								  minX,maxX,minY,maxY,minZ,maxZ, 
								  partIDInSubVolume,
								  parcelIDInSubVolume 	  		);										
					
								  parcelParticleListList[npartI] = parcelParticleList;
								  
								  					
					//Info << " parcelParticleListList " << parcelParticleListList[npartI] << endl;					
					
					/*
					if ( npartI == 0 )
					Foam::createParcels(mesh,particleCloud,parcelCloud,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										
					
										parcelParticleListList[npartI] = parcelParticleList;
										
					if ( npartI == 1 )
					Foam::createParcels(mesh,particleCloud,parcelCloud_2,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										

										parcelParticleListList[npartI] = parcelParticleList;
					
					if ( npartI == 2 )
					Foam::createParcels(mesh,particleCloud,parcelCloud_3,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										

										parcelParticleListList[npartI] = parcelParticleList;

					if ( npartI == 3 )
					Foam::createParcels(mesh,particleCloud,parcelCloud_4,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										

										parcelParticleListList[npartI] = parcelParticleList;
										
					if ( npartI == 4 )
					Foam::createParcels(mesh,particleCloud,parcelCloud_5,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										

										parcelParticleListList[npartI] = parcelParticleList;		
																		
					//Info << " parcelParticleListList " << parcelParticleListList[npartI] << endl;	
					*/
										
				}
			}							
															
			for(int FilterWidth = minFilterWidth; FilterWidth <= maxFilterWidth; FilterWidth += FilterIncrement )
			{
				Info << " " << endl;
				int Filter = 2*(FilterWidth+1)+1;
				Info << "Filter size = " << Filter << "X" << Filter << "X" << Filter << endl;

				if(!Pstream::parRun())
				{					
					filteringEulerianVariables
					(
						args, 
						runTime, 
						mesh,
						StencilListFilter,
						FilterWidth,
						voidfraction,
						U,
						p,
						baralpf,				
						tildeUf,
						barPg,
						Up,
						tildeUp,
						EulerEulerFiltering, 
						EulerianVslipBin
					);
				}else
				{
					Info << " WARNING:: Filtering subroutine is not working in parallel computation" << endl;
					readFilteredEulerianVariables
					(
						args, 
						runTime, 
						mesh,
						baralpf,
						tildeUf,
						barPg,
						FilterWidth
					);
				}
				
				Info << " " << endl;				
				Pout << "Fluid coarsening..." << endl;			
				coarsening(particleCloud,
						   rho,
						   voidfraction,
						   U,
						   tildeUf,
						   p,
						   barPg,
						   gii,
						   particleList,
				   		   verbose,
						   partIDInSubVolume	); // False for parcel	
												
				if(EulerEulerFiltering)
				{
					
					filteredEulerEulerDragCoefficient
					(
						particleCloud,
						verbose,
						mesh,
						StencilListFilter,
						FilterWidth,					
						voidfraction,
						U,
						Up,						
						baralpf,
						tildeUf,
						tildeUp,					
						barPg,
						p,
						rho,
						EulerianDragForce,
						ResolvedEulerianDragForce,
						Eulerian_gii
					);
					
				}
								
				if(particleCoarsening)
				{
					//forAll(npartParticleCoarsening,npartI)
					for ( int npartI = 0; npartI < npartParticleCoarsening.size(); npartI++ )
					{
						// Particle coarsening active or not														
							/*
							ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud,
								verboseParticleCoarsening,voidfraction,rho,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting,
								partIDInSubVolume,
								parcelIDInSubVolume	);
							*/
						
						Info << " " << endl;
						Pout << "Parcel coarsening..." << endl;
						 coarsening(particleCloud,
									rho,
									voidfraction,
									U,
									tildeUf,
									p,
									barPg,
									parcelgii,
									parcelParticleListList[npartI],
									verbose,
									parcelIDInSubVolume				); 
																
																							
					      applyBins(parcelCloud,
									mesh,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									voidfraction,baralpf,
									U,tildeUf,
									Up,tildeUp,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI],
									false,							 // Lagrangian filtering --> false, Eulerian filtering --> true
									EulerianVslipBin,
									parcelIDInSubVolume					); 
												
						/*
						if ( npartI == 0 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
																
							Foam::applyBins(parcelCloud,
									mesh,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );
									

																
						}	
						
						if ( npartI == 1 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud_2,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
								
							Foam::applyBins(parcelCloud_2,
									mesh,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );						
									
						}
						
						if ( npartI == 2 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud_3,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
								
							Foam::applyBins(parcelCloud_3,
									mesh,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );

						}
						
						if ( npartI == 3 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud_4,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
								
							Foam::applyBins(parcelCloud_4,
									mesh,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );

						}
						
						if ( npartI == 4 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud_5,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
								
							Foam::applyBins(parcelCloud_5,
									mesh,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );								

								     
						}
						*/
						
					}
				}							

				Foam::applyBins
				(				
					particleCloud,
					mesh,
					FilterWidth,
					maxFilterWidth,
					nBin,
					minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
					rhop,rho,
					voidfraction,baralpf,
					U,tildeUf,
					Up,tildeUp,
					gii,
					giiCondalppFilter,					
					velSlipCondalppFilter,
					taupVelSlipCondalppFilter,
					numberCondalppFilter,
					VelSlipJointgiiCondalppFilter,
					TaupVelSlipJointgiiCondalppFilter,
					numberVelSlipJointgiiCondalppFilter,
					numberTaupVelSlipJointgiiCondalppFilter,
					false,
					EulerianVslipBin,
					partIDInSubVolume				
				);
				
				if(EulerEulerFiltering)
				{
					Foam::applyBins
					(				
						particleCloud,
						mesh,
						FilterWidth,
						maxFilterWidth,
						nBin,
						minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
						rhop,rho,
						voidfraction,baralpf,
						U,tildeUf,
						Up,tildeUp,
						Eulerian_gii,
						Eulerian_giiCondalppFilter,					
						Eulerian_velSlipCondalppFilter,
						Eulerian_taupVelSlipCondalppFilter,
						Eulerian_numberCondalppFilter,
						Eulerian_VelSlipJointgiiCondalppFilter,
						Eulerian_TaupVelSlipJointgiiCondalppFilter,
						Eulerian_numberVelSlipJointgiiCondalppFilter,
						Eulerian_numberTaupVelSlipJointgiiCondalppFilter,
						true,
						false, // EulerianVslipBin --> false				
						partIDInSubVolume
					);					
				}						
			
			}	
		
		}	

	}
	
	//if(Pstream::master())
	//{	
		// Write bins		
		writeBins
		(
		        runTime,
			mesh,
			nBin,
			minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,		
			outputRelativePath,
			minFilterWidth, 
			maxFilterWidth,
			FilterIncrement,
			-1,				// Just for fluid coarsening results
			giiCondalppFilter,
			velSlipCondalppFilter,
			taupVelSlipCondalppFilter,
			numberCondalppFilter,
			VelSlipJointgiiCondalppFilter,
			TaupVelSlipJointgiiCondalppFilter,
			numberVelSlipJointgiiCondalppFilter,
			numberTaupVelSlipJointgiiCondalppFilter,
			verbose,
			false,  // Lagrangian filtering --> false, Eulerian filtering --> true						
		        EulerianVslipBin 
		);
		
		if(particleCoarsening)
		{
			//forAll(npartParticleCoarsening,npartI)
			for ( int npartI = 0; npartI < npartParticleCoarsening.size(); npartI++ )		
			{
				writeBins
				(
					runTime,
	  				mesh,
					nBin,
					minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,		
					outputRelativePath,
					minFilterWidth, 
					maxFilterWidth,
					FilterIncrement,
					npartParticleCoarsening[npartI],
					NparcelgiiCondalppFilter[npartI],
					NparcelVelSlipCondalppFilter[npartI],
					NparcelTaupVelSlipCondalppFilter[npartI],
					NparcelnumberCondalppFilter[npartI],
					NparcelVelSlipJointgiiCondalppFilter[npartI],
					NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
					NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
					NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI],				
					verboseParticleCoarsening,
					false,  // Lagrangian filtering --> false, Eulerian filtering --> true							
				        EulerianVslipBin
				);			
			}
		}					

		if(EulerEulerFiltering)
		{
			Foam::writeBins
			(
				runTime,
				mesh,
				nBin,
				minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,		
				outputRelativePath,
				minFilterWidth, 
				maxFilterWidth,
				FilterIncrement,
				-1,				// Just for fluid coarsening results
				Eulerian_giiCondalppFilter,
				Eulerian_velSlipCondalppFilter,
				Eulerian_taupVelSlipCondalppFilter,
				Eulerian_numberCondalppFilter,
				Eulerian_VelSlipJointgiiCondalppFilter,
				Eulerian_TaupVelSlipJointgiiCondalppFilter,
				Eulerian_numberVelSlipJointgiiCondalppFilter,
				Eulerian_numberTaupVelSlipJointgiiCondalppFilter,
				verbose,
				true,  // Lagrangian filtering --> false, Eulerian filtering --> true							
			        false 
			);	
		}
	//}	
		
	particleCloud.dataExchangeM().destroy(positions,3);
	particleCloud.dataExchangeM().destroy(velocities,3);

	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
		
    return 0;
}


// ************************************************************************* //
