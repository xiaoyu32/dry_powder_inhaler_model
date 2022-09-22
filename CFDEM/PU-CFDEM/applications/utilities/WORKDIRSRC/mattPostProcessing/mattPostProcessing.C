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
	
	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
					
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
	
	// Open file for results
	OFstream* sPtrVelSlip;
	sPtrVelSlip =  new OFstream("velSlip");
    	*sPtrVelSlip  << "#Time		< <u_f> - v_p >_p	< u_f > - < v_p >_p" << endl; 			

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

		// Read gas voidfraction
		IOobject alpfheader
		(
		   "phiFluid",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		   
		Info<< " 	Reading alpf" << endl;
		volScalarField alpf(alpfheader,mesh);
		
		// Domain average gas velocity
		vector domainAveUf(0,0,0);
		domainAveUf = fvc::domainIntegrate(alpf*U).value()/fvc::domainIntegrate(alpf).value(); 

		// Drift velocity
                vector domainAveDrift(0,0,0);
                domainAveDrift = fvc::domainIntegrate((1.-alpf)*U).value()/fvc::domainIntegrate(1.-alpf).value() - domainAveUf;

                Info << " Drift velocity = " << domainAveDrift << endl;
		
		// Create gas velocity interpolation variable
		interpolationCellPoint<vector> Uf_xp(U);					

		// < u_f|x_p - v_p >
		vector velSlip(0,0,0);
		// < u_f - v_p >
		vector velSlipNonInt(0,0,0);
		// < < u_f > - v_p >
		vector tildeVelSlip(0,0,0);
	
		// Number or particles
		int nP(0);
		int just_integer(-1);
		int countCellID(0);
		
		if ( runTime.timeName() != "0" )
		{

			// Read particle positions 
			fileName H(runTime.timeName()/"positions");
			Info << "" << endl;
			Info << " Opening file " << H << endl;
	        	// set file pointer
	        	string HH=string(H);
	        	const char * particleFilePath=HH.c_str();
	        	ifstream* inputPtr;
	        	inputPtr = new ifstream(particleFilePath);

        		// read data
        		string just_read;
			int countLine(0);
			
			std::getline(*inputPtr,just_read);
			while(countLine<98)
			{
				*inputPtr >>just_read;
				//Info << " just_read " << just_read << endl;
				countLine++;
			}

			*inputPtr >> nP;		// Read number of particle
			Info << " Number of particles = " << nP << endl;			
			
			*inputPtr >> just_read;		// Read empty paranthesis
        		
			for(int index = 0;index < nP; ++index)
        		{
        		    *inputPtr >> just_read >> positions[index][0] >> positions[index][1] >> positions[index][2] >>  just_read >> just_integer ;
        		    // Particle diameter
			    radii[index][0] = 145.e-06/2.;
			    //Info << positions[index][0] << positions[index][1] << positions[index][2] << endl;
			}

			// clean up inputStream
    			delete inputPtr;

			// Read particle velocities
			fileName HVel(runTime.timeName()/"Up");
			Info << "" << endl;
			Info << " Opening file " << HVel << endl;
	        	// set file pointer
	        	string HHVel=string(HVel);
	        	const char * particleFilePathVel=HHVel.c_str();
	        	ifstream* inputPtrVel;
	        	inputPtrVel = new ifstream(particleFilePathVel);

        		// read data
			countLine=0;
			std::getline(*inputPtrVel,just_read);
			while(countLine<98)
			{
				*inputPtrVel >>just_read;
				//Info << " just_read " << just_read << endl;
				countLine++;
			}
			*inputPtrVel >> nP;		// Read number of particle
			Info << " Number of particles = " << nP << endl;			
			*inputPtrVel >> just_read;		// Read empty paranthesis
        		
			for(int index = 0;index < nP; ++index)
        		{
        		    *inputPtrVel >> just_read >> velocities[index][0] >> velocities[index][1] >> velocities[index][2] >>  just_read ;
        		    //Info << velocities[index][0] << velocities[index][1] << velocities[index][2] << endl;
			}

			// clean up inputStream
    			delete inputPtrVel;

			particleCloud.setRadii(radii);
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
					if(particleCloud.particleCell(index)>-1) Info << " Uf     = " << U[particleCloud.particleCell(index)] << endl;
					if(particleCloud.particleCell(index)>-1) Info << " Uf@Xp  = " << Uf_xp.interpolate(particleCloud.position(index),particleCloud.particleCell(index));
					Info << "" << endl;
				}	
			}
			

			label cellID(-1);
			countCellID = 0;
			for( int partI = 0; partI < particleCloud.numberOfParticles(); partI++ )
			{
			     tildeVelSlip  += particleCloud.velocity(partI);	
			     cellID = particleCloud.particleCell(partI);
			     if(cellID>-1)
			     {	
				countCellID++;
				velSlip       += domainAveUf - particleCloud.velocity(partI); //Uf_xp.interpolate(particleCloud.position(partI),cellID) - particleCloud.velocity(partI);
				velSlipNonInt += U[cellID] - particleCloud.velocity(partI);					
			     }
			} 
						
		}
		
		Info << "" <<endl;
		Info << " Number of mapped particle = " << countCellID  << endl;
		Info << " < <u_f> - v_p >_p = " << velSlip/nP	 	<< endl;
		Info << " < u_f > - < v_p >_p = " << domainAveUf - tildeVelSlip/nP	<< endl;			
		Info << "" <<endl;
		
    		*sPtrVelSlip 	<< U.mesh().time().value()      << tab 
                		<< velSlip[0]/nP 	<< tab << velSlip[1]/nP 	<< tab << velSlip[2]/nP 	<< tab
				<< domainAveUf[0] - tildeVelSlip[0]/nP 	<< tab << domainAveUf[1] - tildeVelSlip[1]/nP 	<< tab << domainAveUf[2] - tildeVelSlip[2]/nP	<< endl;	
		
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
