/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPiso

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "cfdemCloud.H"
#include "cfdemCloudCharge.H"
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
    CfdemCloudCharge particleCloud(mesh);

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
        particleCloud.evolveCharge(Ef,rhoe,voidfraction,Us,U);
        //particleCloud.evolve(voidfraction,Us,U);
	
        Info << "update Ksl.internalField()" << endl;
        Ksl.oldTime().internalField() = particleCloud.momCoupleM(0).impMomSource();
        particleCloud.smoothingM().smoothen(Ksl);
        Ksl.correctBoundaryConditions();

        // Update continuous phase volume fraction field
        voidfraction.correctBoundaryConditions();
        surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
          
        // Update bcs for Us
        Us.correctBoundaryConditions();
 
	// Source term
	fsource.internalField() = particleCloud.momCoupleM(1).expMomSource();
	fsource.correctBoundaryConditions();
	
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
	
	//solve electrical field
	#include "phiEqn.H"
	particleCloud.evolveElectricField( Ef );

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

        Info << "Getting collisional stress data..." << endl;
        if (runTime.outputTime() && calcStress) 
	{
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


