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

#include "fluidParticleCoarsening.H"
#include "domainFilter.H"
#include "readVariables.H"
#include "dragForce.H"
#include "fluidCoarsening.H"
#include "createParcel.H"
#include "coarseningPar.H"
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
	//particleCloud.reAllocArrays();
	particleCloud.reAllocArraysPost();
    
    	// Create parcel cloud
    	cfdemCloud parcelCloud(mesh);
    	//parcelCloud.reAllocArrays();
	parcelCloud.reAllocArraysPost();

	double **positions;
	double **velocities;
	double **radii;
	double **cellID;

	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
    	//particleCloud.dataExchangeM().allocateArray(radii,0.,1);
    	//particleCloud.dataExchangeM().allocateArray(cellID,0.,30);
  
    	particleCloud.get_radii(radii);
	particleCloud.get_cellIDs(cellID);
    
	// Parcel variables
	double **parcelPositions;
        double **parcelVelocities;
        double **parcelRadii;
        double **parcelCellID;

        parcelCloud.dataExchangeM().allocateArray(parcelPositions,0.,3);
        parcelCloud.dataExchangeM().allocateArray(parcelVelocities,0.,3);
        //particleCloud.dataExchangeM().allocateArray(radii,0.,1);
        //particleCloud.dataExchangeM().allocateArray(cellID,0.,30);

        parcelCloud.get_radii(parcelRadii);
        parcelCloud.get_cellIDs(parcelCellID);
	
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
	int nBin(10);
	scalar maxalpp(0);
	scalar minalpp(0.64);
	// Binning with VelSlip, taup |Vr|
	scalar minVelSlip(-1);
	scalar maxVelSlip( 1);
	scalar minTaupVelSlip(0);
	scalar maxTaupVelSlip(0.5);
	//
    if(propsDict.found("nBin"))             nBin = readScalar(propsDict.lookup("nBin"));
    if(propsDict.found("Max_alpp"))         maxalpp = readScalar(propsDict.lookup("Max_alpp"));
    if(propsDict.found("Min_alpp"))         minalpp = readScalar(propsDict.lookup("Min_alpp"));
	if(propsDict.found("Min_VelSlip"))      minVelSlip = readScalar(propsDict.lookup("Min_VelSlip"));
	if(propsDict.found("Max_VelSlip"))      maxVelSlip = readScalar(propsDict.lookup("Max_VelSlip"));
	if(propsDict.found("Min_TaupVelSlip"))  minTaupVelSlip = readScalar(propsDict.lookup("Min_TaupVelSlip"));
	if(propsDict.found("Max_TaupVelSlip"))  maxTaupVelSlip = readScalar(propsDict.lookup("Max_TaupVelSlip"));
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
	
	// Data points
	// ANNpointArray dataPts;
	// Query points
	// ANNpoint queryPt;
	// Search structure
	// ANNkd_tree* kdTree;
	
	// ######################### NEW VARIABLES (Parallel compuatatio) #######################
					
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Info << " " << endl;
		Info << "Time = " << runTime.timeName() << endl;

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
		
            //particleCloud.get_cellIDs(cellID);
            //particleCloud.setRadii(radii);
	    //particleCloud.get_radii(radii);
            particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
	    particleCloud.setPos(positions);
	    particleCloud.setVel(velocities);
            
            // Create particle list
            labelList particleL(1,(0));
            labelListList particleList(particleCloud.numberOfParticles(),particleL);
            
            // Nunber of particles in sub-volume
            int nPartInSubVolume = 0;
            scalar totalPartVol = 0;
						
			// Particle IDs in sub-volume
			for(int index = 0; index < particleCloud.numberOfParticles(); index++)
	        {
				label cellI = particleCloud.cellIDs()[index][0];					
				if(cellI > -1)
				{
					totalPartVol += 4./3. * constant::mathematical::pi
							      * particleCloud.radius(index) 
							      * particleCloud.radius(index)
							      * particleCloud.radius(index); 
							      
					// Pout << ii <<  " cellI " << cellI  << " partIDInSubVolume[ii] "	<< partIDInSubVolume[ii] << endl;
                    particleList[nPartInSubVolume][0] = index;
                    nPartInSubVolume++;
				}
			}	
					
			// Resize partID list in sub-volume
            particleList.resize(nPartInSubVolume);
			alppInSubVolume = totalPartVol / domainVol;
			if(Pstream::parRun()) Pout << " Number of particles in sub-volume = " << nPartInSubVolume <<endl;
			if(Pstream::parRun()) Pout << " Solid volume fraction in sub-volume = " << alppInSubVolume <<endl;

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
				particleList
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
					particleList			
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
            // Parcel List
			labelListList parcelList;
            
			if(particleCoarsening)
			{
				for ( int npartI = 0; npartI < npartParticleCoarsening.size(); npartI++ )
				{
					Info << " " << endl;
					Info << " Particle coarsening: "<< endl;
					
					// Number of parcel
                    int nParcel =  particleList.size();
                    labelList particleIDList(npartParticleCoarsening[npartI],(0));
					labelListList parcelParticleList(nParcel,particleIDList);
					parcelParticleListList.setSize(npartParticleCoarsening.size(),parcelParticleList);
					
                    labelList parcelL(1,(0));
                    //parcelList.setSize(nParcel,parcelL);
                    // Find cellID algorithm moves some particle from sub-volume to other; Just for security
                    parcelList.setSize(nParcel+1000,parcelL);
                    // parcelList.setSize(nParcel,parcelL);

					Info << " " << endl;				
					Pout << "Creating parcels..." << endl;
					createParcels(
                                  mesh,
                                  particleCloud,
                                  parcelCloud,
                                  npartParticleCoarsening[npartI],
                                  verboseParticleCoarsening,
                                  minX,
                                  maxX,
                                  minY,
                                  maxY,
                                  minZ,
                                  maxZ,
                                  // dataPts,
                                  // queryPt,
                                  // kdTree,
				  parcelPositions,
        			  parcelVelocities,
        			  parcelRadii,
        			  parcelCellID,
                                  particleList,
                                  parcelList,
                                  parcelParticleList);
                    
                                  Pout << "Difference between parcel list and parcelID " << parcelList.size() - parcelParticleList.size() << endl;
					
								  parcelParticleListList[npartI] = parcelParticleList;
                    
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
				coarseningPar(particleCloud,
                              rho,
                              voidfraction,
                              U,
                              tildeUf,
                              p,
                              barPg,
                              gii,
                              particleList,
                              verbose       );
												
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
					for ( int npartI = 0; npartI < npartParticleCoarsening.size(); npartI++ )
					{
						
						Info << " " << endl;
						Pout << "Parcel coarsening..." << endl;
						 coarseningPar(particleCloud,
                                       rho,
                                       voidfraction,
                                       U,
                                       tildeUf,
                                       p,
                                       barPg,
                                       parcelgii,
                                       parcelParticleListList[npartI],
                                       verbose);
																                        
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
                                    parcelList);
						
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
                    particleList
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
                        particleList
                     );
				}						
			
			}	
		
		}	

	}
	
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
