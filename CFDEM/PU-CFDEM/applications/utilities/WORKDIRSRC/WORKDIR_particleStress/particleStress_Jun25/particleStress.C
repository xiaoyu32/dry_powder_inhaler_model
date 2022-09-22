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

#include "calcCollisionalForce.H"

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

	// Create particle cloud			
	cfdemCloud particleCloud(mesh);
	
	particleCloud.reAllocArrays();

	double **positions;
	double **velocities;
	double **omegas;
	double **radii;
	double **cellID;
	
	double **forces;

	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.dataExchangeM().allocateArray(omegas,0.,3);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
	
	particleCloud.dataExchangeM().allocateArray(forces,0.,3);
				
	// Read dictionary + dictionaryProps 
	const dictionary dict(particleCloud.couplingProperties());
	int maxNumberOfParticles(0);
	const dictionary oneWayVTKPropsDict(dict.subDict("oneWayVTKProps"));
	maxNumberOfParticles=readScalar(oneWayVTKPropsDict.lookup("maxNumberOfParticles"));
	// Read particleStress Sub-dictionary
	const dictionary ParticleStressPropsDict(dict.subDict("ParticleStressProps"));
	const scalar rhop(readScalar(ParticleStressPropsDict.lookup("rhoParticle")));					
	bool verbose(false);
	if(ParticleStressPropsDict.found("verbose")) verbose = true;	
	
	const scalar collisionDp(readScalar(ParticleStressPropsDict.lookup("collisionDp")));

	// DEM input file parameters
	fileName DEMinputFileRelativePath("../DEM");
	fileName DEMinputFilename(ParticleStressPropsDict.lookup("DEMinputFilename"));

	// Check DEM input file exist or not
	if( !isFile(mesh.time().path()/DEMinputFileRelativePath/DEMinputFilename) )
	{
		FatalError << " DEM input file does not exist " << abort(FatalError);
	}
	else
	{
		Info << " " << endl;
		Info << "DEM input file is " << mesh.time().path()/DEMinputFileRelativePath/DEMinputFilename << endl;
	}

	// Collision parameters
	#include "collision_parameters.H"
	
	// Total volume of domain
	scalar domainVol(0);
	forAll(mesh.C(),cellI)
	{
		domainVol +=mesh.V()[cellI];
	}
	Info << "Domain volume[m^3] = " << domainVol << endl;	
								
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();

		if ( runTime.timeName() != "0" )
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

			particleCloud.dataExchangeM().getData("f","vector-atom",forces,count);
			Info<< " 	Reading forces on particles" << endl;
			Info<< " " << endl;				
		
			particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
			particleCloud.setPos(positions);
			particleCloud.setVel(velocities);
			particleCloud.setOmega(omegas);
			
			if(verbose)
			{
				int index = 0;
				Info << "" << endl;
				Info << " index  = " << index << endl;
				Info << " rp     = " << particleCloud.radius(index) << endl;
				Info << " Vp     = " << particleCloud.velocity(index) << endl;
				Info << " Omegap = " << particleCloud.omega(index) << endl;				
				Info << " Xp     = " << particleCloud.position(index) << endl;
				Info << " CellID = " << particleCloud.particleCell(index) << endl;
				Info << " Fp     = " << forces[index][0] << " " << forces[index][1] << " " << forces[index][2] << endl;
				Info << "" << endl;
			}				
			
			// Neighboring list parameters
			int k = 20; 
			if(particleCloud.numberOfParticles()<k) k = particleCloud.numberOfParticles();
			
			// Dimensions, exact OR approximate  
			int dim=3; double eps = 0;

			// Number of points
			int nPts;
			nPts =  particleCloud.numberOfParticles();		       		

			// Data points
			ANNpointArray dataPts;
			// Query points
			ANNpoint queryPt;

			ANNidxArray	nnIdx;			// 	near neighbour indices
			ANNdistArray 	dists;			//	near neighbour distances
			ANNkd_tree* 	kdTree;			//	search structure

			// Allocate 
			queryPt = annAllocPt(dim);
			dataPts = annAllocPts(nPts, dim);
			nnIdx = new ANNidx[k];
			dists = new ANNdist[k];

			// Particle kinetic stresses
			symmTensor sigma_kin(0,0,0,0,0,0);
			
			// Particle collisional stresses
			symmTensor sigma_coll(0,0,0,0,0,0);
			
			// Mean particle velocity
			scalar meanVpx(0);scalar meanVpy(0);scalar meanVpz(0); 
			
			// Create particle list
			for(int index = 0; index <  particleCloud.numberOfParticles(); index++)
	        	{			
				dataPts[index][0] = particleCloud.position(index).x();		
				dataPts[index][1] = particleCloud.position(index).y();		
				dataPts[index][2] = particleCloud.position(index).z();
				
				// Summation to calculate mean particle velocity
				meanVpx +=particleCloud.velocity(index).x();
				meanVpy +=particleCloud.velocity(index).y();
				meanVpz +=particleCloud.velocity(index).z();

				for (int dir=0;dir<3;dir++)
				{	
					fcoll[index][dir] = 0;
					 ftan[index][dir] = 0;		
 					 fcap[index][dir] = 0;
 					fvisc[index][dir] = 0;
				}

			}
			
			// Mean particle velocity
			meanVpx /= particleCloud.numberOfParticles();
			meanVpy /= particleCloud.numberOfParticles();
			meanVpz /= particleCloud.numberOfParticles();	
			
			Info << " " << endl;
			Info << " Mean particle velocities " << endl;
			Info << " <u_p,x> = " << meanVpx << endl;
			Info << " <u_p,y> = " << meanVpy << endl;
			Info << " <u_p,z> = " << meanVpz << endl;					
			Info << " " << endl;

			// Create Tree structure
			kdTree = new ANNkd_tree(dataPts, nPts, dim);

			vector pos_i(0,0,0);
			vector pos_j(0,0,0);
			scalar dist(0);
			
			// Particle mass
			scalar mass_p(0);
					
			for(int index = 0; index <  particleCloud.numberOfParticles(); index++)		
			{
				// Info << " Omegap = " << particleCloud.omega(index) << endl;	// To test 
				
				scalar sqRad = collisionDp*particleCloud.radius(index);
				
				queryPt[0] = particleCloud.position(index).x();
				queryPt[1] = particleCloud.position(index).y();
				queryPt[2] = particleCloud.position(index).z();

				kdTree->annkFRSearch(
			                		queryPt,			// query point					
							sqRad,				// squared radius
							k,				// number of the near neighbours to return
							nnIdx,				// nearest neighbor array
							dists,				// dist to near neighbours
							eps			);					

				for (int index_j = 1; index_j < k; index_j++) // k=0 particle itself
				{
					// Calculate collision forces
					calcForce(			  
							particleCloud,
						      collisionModelI, 
								index,
						       nnIdx[index_j],     
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
							   sigma_coll  	);				  						
					
				}
				
				// Particle kinetic stress
				mass_p = 4./3.*rhop*constant::mathematical::pi*particleCloud.radius(index)*particleCloud.radius(index)*particleCloud.radius(index);	
				sigma_kin.xx() += mass_p * ( particleCloud.velocity(index).x() - meanVpx ) * ( particleCloud.velocity(index).x() - meanVpx ); 
				sigma_kin.xy() += mass_p * ( particleCloud.velocity(index).x() - meanVpx ) * ( particleCloud.velocity(index).y() - meanVpy ); 
				sigma_kin.xz() += mass_p * ( particleCloud.velocity(index).x() - meanVpx ) * ( particleCloud.velocity(index).z() - meanVpz ); 
				sigma_kin.yy() += mass_p * ( particleCloud.velocity(index).y() - meanVpy ) * ( particleCloud.velocity(index).y() - meanVpy ); 
				sigma_kin.yz() += mass_p * ( particleCloud.velocity(index).y() - meanVpy ) * ( particleCloud.velocity(index).z() - meanVpz ); 	
				sigma_kin.zz() += mass_p * ( particleCloud.velocity(index).z() - meanVpz ) * ( particleCloud.velocity(index).z() - meanVpz ); 

			}
						
			if(verbose)
			{
				int index = 0;
				Info << "" << endl;
				Info << " Ftot = " << fcoll[index][0] << " " << fcoll[index][1] << " " << fcoll[index][2] << endl;
				if(cohesion) Info << " Fcoh = " <<  fcoh[index][0] << " " <<  fcoh[index][1] << " " <<  fcoh[index][2] << endl;
				Info << "" << endl;
			}
			
			Info << " Particle kinetic stresses" << endl;
			Info << " sigma_kin_xx= " << sigma_kin.xx()/domainVol << endl;
			Info << " sigma_kin_yy= " << sigma_kin.yy()/domainVol << endl;
			Info << " sigma_kin_zz= " << sigma_kin.zz()/domainVol << endl;			
			Info << " sigma_kin_xy= " << sigma_kin.xy()/domainVol << endl;
			Info << " sigma_kin_xz= " << sigma_kin.xz()/domainVol << endl;
			Info << " sigma_kin_yz= " << sigma_kin.yz()/domainVol << endl;
			
			Info << "" << endl;
			/*
			Info << " Particle collisional stresses" << endl;
			Info << " sigma_coll_xx= " << sigma_coll.xx()/2./domainVol << endl;
			Info << " sigma_coll_xy= " << sigma_coll.xy()/2./domainVol << endl;
			Info << " sigma_coll_xz= " << sigma_coll.xz()/2./domainVol << endl;
			Info << " sigma_coll_yy= " << sigma_coll.yy()/2./domainVol << endl;
			Info << " sigma_coll_yz= " << sigma_coll.yz()/2./domainVol << endl;
			Info << " sigma_coll_zz= " << sigma_coll.zz()/2./domainVol << endl;		
			*/
			Info << " Particle collisional stresses" << endl;
			Info << " sigma_coll_xx= " << sigma_coll.xx()/domainVol << endl;
			Info << " sigma_coll_yy= " << sigma_coll.yy()/domainVol << endl;
			Info << " sigma_coll_zz= " << sigma_coll.zz()/domainVol << endl;						
			Info << " sigma_coll_xy= " << sigma_coll.xy()/domainVol << endl;
			Info << " sigma_coll_xz= " << sigma_coll.xz()/domainVol << endl;			
			Info << " sigma_coll_yz= " << sigma_coll.yz()/domainVol << endl;
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
