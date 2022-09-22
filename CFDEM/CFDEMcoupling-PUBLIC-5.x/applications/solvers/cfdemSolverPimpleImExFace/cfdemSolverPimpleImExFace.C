/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
                                Copyright (C) 2013-     Graz University of
                                                        Technology, IPPT
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
    cfdemSolverPimpleImExFace

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pimpleFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.

    This solver can handle MRFZones, as well as explicitPorousZones via the fvOptions
    feature introduced in OF 2.2.1

    This solver can handle cylPorousZones needed by some solvers
    to model rotating porous regions

    It calculates the pressure and momentum exchange terms on the faces.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"

#include "OFversion.H"
#if defined(version30)
    #include "turbulentTransportModel.H"
    #include "pisoControl.H"
#else
    #include "turbulenceModel.H"
#endif

#if defined(versionv1606plus) || defined(version40)
    #include "fvOptions.H"
#else
    #include "fvIOoptionList.H"
#endif

#include "fixedFluxPressureFvPatchScalarField.H"

#include "cfdemCloud.H"
#if defined(biDisperseModels)
    #include "cfdemCloudBiDisperse.H"
    #include "cfdemCloudBiDisperseRotation.H"
#endif
#include "implicitCouple.H"
#include "forceModel.H"
#include "clockModel.H"
#include "smoothingModel.H"

#include "pimpleControl.H"
#include "interpolationCellPoint.H"

#include "scalarTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #if defined(version30)
        #include "createTimeControls.H"
        #include "createMRF.H"
    #endif
    #include "createFields.H"

    //fvOptions and other zones for adding sources
    #include "createFvOptions.H"
    #if defined(biDisperseModels)
    #include "createFieldsBiDisperse.H"
    #endif
    #include "createFieldsAddOn.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    #if defined(biDisperseModels)
    cfdemCloudBiDisperseRotation    particleCloud(mesh);
    #else
    cfdemCloud particleCloud(mesh);
    #endif
    IOdictionary couplingProperties = particleCloud.couplingProperties();
    bool writeOnlyParticles=false;
    if (couplingProperties.found("writeOnlyParticles") && !particleCloud.solveFlow())
    {
        writeOnlyParticles=true;
        Info << "WARNING: will only write particle data, no fluid data! " << endl;
    }

    #if defined(biDisperseModels)
    bool doBiDisperseUpdates = true;
    if (couplingProperties.found("skipBiDisperseUpdates"))
    {
        doBiDisperseUpdates=false;
        Info << "WARNING: will not update Euler fields necessary for bi-disperse drag models! " << endl;
    }
    #endif

    #include "checkModelType.H"
    #include "checkImExCoupleM.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    autoPtr<scalarTransportModel> stm
    (
        scalarTransportModel::New(particleCloud.couplingProperties(),particleCloud)
    );
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
        particleCloud.clockM().start(1,"Global");


        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #if defined(version30)
            #include "readTimeControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"
        #else
            #include "readPISOControls.H"
            #include "CourantNo.H"
        #endif

        // do particle stuff
        #include "evolveCloud.H"

        //Scalar transport if desired. Use "none" (noTransport) if no scalar transport is desired
        stm().update();

        // Pressure-velocity PIMPLE Algorithm
        if(particleCloud.solveFlow())
        {
            while (pimple.loop())
            {

                if(ignoreFinalPimpleIteration)   mesh.data::remove("finalIteration"); //remove the finalIteration
                                                                                                                //keyword if phyiscal timestepping
                                                                                                                //not desired

                #include "UEqn.H"

                #if defined(version30)
                    while (pimple.correct())
                #else
                    for (int corr=0; corr<nCorr; corr++)
                #endif
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    laminarTransport.correct();
                    turbulence->correct();
                }
            }
        }

        //update density
        fvOptions.correct(rho);

        #if defined(biDisperseModels)
            if(!writeOnlyParticles) //only write Euler fields if necessary
            {
                if(runTime.write() && doBiDisperseUpdates)
                {
                    uP1.write();
                    uP2.write();
                    phiP1.write();
                    phiP2.write();
                    dSauter.write();
                }
            }
        #else
            runTime.write();
        #endif

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
