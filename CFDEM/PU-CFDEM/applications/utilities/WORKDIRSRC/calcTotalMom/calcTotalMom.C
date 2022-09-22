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

#include "calcTotalMom.H"
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
    	
    double **velocities;
    double **radii;

    particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
    particleCloud.get_radii(radii); 
    
    forAll(timeDirs, timeI)
    {

	runTime.setTime(timeDirs[timeI], timeI);
	if ( runTime.timeName() != "0" )
	{
        	
		Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;

        	mesh.readUpdate();

		// Read gas density
		IOobject rhoheader
		(
		   "rho",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		volScalarField rho_f(rhoheader,mesh);	

		// Read volume fraction of gas
		IOobject voidfractionheader
		(
		   "voidfraction",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);

		volScalarField alp_f(voidfractionheader,mesh);	

		// Read gas velocity 
		IOobject Uheader
		(
		   "U",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);

		volVectorField U_f(Uheader,mesh);

                
        	int count = runTime.value() / particleCloud.dataExchangeM().DEMts();				
		Info<< " " << endl;
		particleCloud.dataExchangeM().getData("v","vector-atom",velocities,count);
		Info<< " 	Reading particle velocities" << endl;
		Info<< " " << endl;
		
		particleCloud.dataExchangeM().getData("radius","scalar-atom",radii,count);
		Info<< " 	Reading particle radius" << endl;		
		Info<< " " << endl;
		
		vector sumUp(0,0,0);
		scalar ds=2*radii[0][0];
		
		for(int index = 0;index <  particleCloud.numberOfParticles(); index++)
    		{		
			for(int j=0;j<3;j++) sumUp[j] += velocities[index][j];
		}	
		// Calculate vslip
		Info << " Domain-averaged slip velocity = " 
		     << fvc::domainIntegrate(alp_f*U_f).value()/fvc::domainIntegrate(alp_f).value() 
		       -sumUp/particleCloud.numberOfParticles() << endl;

		// Calculate momentum
		Info 	<< " Domain-averaged fluid momentum = " 
	     		<< fvc::domainIntegrate(rho_f*alp_f*U_f).value()
			<< " total particle Vol*Up = " << sumUp*ds*ds*ds*M_PI/6 << endl;	
		
	}	
    }
    
    particleCloud.dataExchangeM().destroy(velocities,3);
	
    return 0;
}


// ************************************************************************* //
