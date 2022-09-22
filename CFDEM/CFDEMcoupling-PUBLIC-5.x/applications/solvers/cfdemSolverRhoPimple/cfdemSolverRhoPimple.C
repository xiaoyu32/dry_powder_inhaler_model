/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2015-2020 Thomas Lichtenegger, JKU Linz, Austria
	Copyright (c) 2020-     Stefan Radl, TU Graz, Austria (modification to work with Sundar's LIGGGHTS and CFDEM version)

Application
    cfdemSolverRhoPimple

Description
    Transient solver for compressible flow using the flexible PIMPLE (PISO-SIMPLE)
    algorithm.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver rhoPimpleFoam in OpenFOAM(R) 4.x,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "cfdemCloudEnergy.H"
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"
#include "thermCondModel.H"
#include "energyModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createFieldsAddOn.H" //creates f and fError field
    #include "createFvOptions.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    cfdemCloudEnergy particleCloud(mesh);
    #include "checkModelType.H"

    turbulence->validate();
    //#include "compressibleCourantNo.H"
    //#include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        particleCloud.clockM().start(1,"Global");
		
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
		
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        // do particle stuff
		
        #include "evolveCloud.H"
		
        //Scalar transport if desired. Use "none" (noTransport) if no scalar transport is desired
        //stm().update(); //XXX to implement if desired. see solver PimpleImExFace
		
#if OPENFOAM_VERSION_MAJOR < 6
        if (pimple.nCorrPIMPLE() <= 1)
#else
        if (pimple.nCorrPimple() <= 1)
#endif
        {
            #include "rhoEqn.H"
        }

        volScalarField rhoeps("rhoeps",rho*voidfraction);
        // Pressure-velocity PIMPLE Algorithm
        if(particleCloud.solveFlow())
        {
			while (pimple.loop())
			{
				#include "UEqn.H"
				#include "EEqn.H"

				// --- Pressure corrector loop
				while (pimple.correct())
				{
					// besides this pEqn, OF offers a "pimple consistent"-option
					#include "pEqn.H"
					rhoeps=rho*voidfraction;
				}

				if (pimple.turbCorr())
				{
					turbulence->correct();
				}
			}
        }
		
        particleCloud.clockM().start(31,"postFlow");
        particleCloud.postFlow();
        particleCloud.clockM().stop("postFlow");

        runTime.write();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
