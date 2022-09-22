/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    DPMFoam

Description
    Transient solver for the coupled transport of a single kinematic particle
    cloud including the effect of the volume fraction of particles on the
    continuous phase.

\*---------------------------------------------------------------------------*/




#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "cfdemCloud.H"
#include "cfdemCloudStress.H"
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"
#include "cyclicFvPatch.H"
#include "pimpleControl.H"

#include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"    
    #include "createFields.H"        
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    //cfdemCloud particleCloud(mesh);
    //cfdemCloudPerio particleCloud(mesh); 
    cfdemCloudStress particleCloud(mesh);

    pimpleControl pimple(mesh);

    volScalarField nufField = particleCloud.turbulence().nu();
    
    Info<< "\nStarting time loop\n" << endl;
     
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "\nTime = " << runTime.timeName() << nl << endl;

        Info << "Evolving particle cloud..." << endl;
        particleCloud.evolve(voidfraction,Us,U);

	Ksl.oldTime().internalField() = particleCloud.momCoupleM(0).impMomSource();
        particleCloud.smoothingM().smoothen(Ksl);
        Ksl.correctBoundaryConditions();

        // Update continuous phase volume fraction field
        voidfraction.correctBoundaryConditions();
        surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);

	// Source term
	fsource.internalField() = particleCloud.momCoupleM(1).expMomSource();
	fsource.correctBoundaryConditions();
	
	Info<<"Cell fsource value:"<<fsource[0][2]<<endl;
	
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }

        // Fully periodic flow
        if(fullyPeriodicFLow)
	{
	     #include "momIO.H"
	} 

        // Mixture velocity correction
        if(mixtureVelocityCorrectionPeriodicFlow)
		particleCloud.calcSendVelocityCorrection
		(
		    rhoParticle,
     		    rho,
    		    voidfraction,
                    U,
     		    Us
		);

        if (runTime.outputTime() && calcStress) 
	{
          Info << "Getting collisional stress data..." << endl;
          Info << "Particles stresses are being written..." << endl;
	  particleCloud.getDEMStressData();
	  particleCloud.writeCollStress(mesh,particleCloud); 	
        }
	
	runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
