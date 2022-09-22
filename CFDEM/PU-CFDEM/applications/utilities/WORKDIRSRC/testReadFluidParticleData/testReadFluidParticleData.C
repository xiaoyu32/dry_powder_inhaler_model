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

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"

#include "cfdemCloud.H"
#include "dataExchangeModel.H"
#include "voidFractionModel.H"
#include "locateModel.H"
#include "averagingModel.H"
#include "momCoupleModel.H"
#include "forceModel.H"
#include "IOModel.H"
#include "interpolationCellPoint.H"

#include "timeSelector.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

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
	double **radii;
	double **cellID;
	
	double **forces;

	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
	
	particleCloud.dataExchangeM().allocateArray(forces,0.,3);
				
	// Read dictionary + dictionaryProps 
	const dictionary dict(particleCloud.couplingProperties());
	const dictionary propsDict(dict.subDict("oneWayVTKProps"));					
	bool verbose = false;
	labelList exIndex(1,0);
	// Particle ID for debuging
	if(propsDict.found("verbose"))
	{
		verbose = true;	
		if(propsDict.found("exIndex")) exIndex = labelList(propsDict.lookup("exIndex"));	
	}	
		
	Info << " exIndex = " << exIndex << endl;
		
	// Read scalar from sub-dict
	scalar rhop(0);
	if(propsDict.found("rhoParticle")) rhop = readScalar(propsDict.lookup("rhoParticle"));
	Info << " Particle density = " << rhop << endl;				

	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();
		
		// Read gas velocity 
		IOobject Uheader
		(
		   "U",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		   
		Info<< " 	Reading U" << endl;
		volVectorField U(Uheader,mesh);
		
		// Create gas velocity interpolation variable
		interpolationCellPoint<vector> Uf_xp(U);

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

			particleCloud.dataExchangeM().getData("f","vector-atom",forces,count);
			Info<< " 	Reading forces on particles" << endl;
			Info<< " " << endl;				
		
			// Locate particles in Eulerian grid
			particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
			particleCloud.setPos(positions);
			particleCloud.setVel(velocities);

			if(verbose)
			{
				//int index = 0;
				forAll(exIndex,ii)
				{
					int index = exIndex[ii];
					
					Info << "" << endl;
					Info << " index  = " << index << endl;
					Info << " rp     = " << particleCloud.radius(index) << endl;
					Info << " Vp     = " << particleCloud.velocity(index) << endl;
					Info << " Xp     = " << particleCloud.position(index) << endl;
					Info << " CellID = " << particleCloud.particleCell(index) << endl;
					Info << " Fp     = " << forces[index][0] << " " << forces[index][1] << " " << forces[index][2] << endl;
					
					int theCellIdOfTheParticle = particleCloud.particleCell(index);
					
					Info << " Uf     = " << U[theCellIdOfTheParticle] << endl;
					Info << " Uf@Xp  = " << Uf_xp.interpolate(particleCloud.position(index),theCellIdOfTheParticle);
					Info << "" << endl;
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
