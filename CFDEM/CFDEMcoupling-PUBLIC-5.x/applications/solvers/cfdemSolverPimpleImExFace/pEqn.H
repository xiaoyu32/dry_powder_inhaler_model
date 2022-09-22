
                volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);
                surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidfraction)); // ChG

                volVectorField HbyA("HbyA", U);

                #if defined(versionv1606plus) || defined(version40)
                    HbyA = rUA*UEqn.H();
                #else
                    HbyA = rUA*UEqn().H();
                #endif

                //if (nCorrSoph <= 1)
                //{
                //    UEqn.clear();
                //}

                phi = voidfractionf*phiByVoidfraction;

                //Compute the provisional flux (NOT superficial, since voidfraction will be multiplied directly in pEqn)
                #ifdef version23
                    surfaceScalarField phiHbyA
                    (
                        "phiHbyA",
                        (
                            (fvc::interpolate(HbyA) & mesh.Sf() )
                            + phicForces                               //explicit contribution
                            + rUAfvoidfraction*fvc::ddtCorr(U, phiByVoidfraction) //correction
                        )
                    );

                    #if defined(version30)
                        MRF.makeRelative(phiByVoidfraction);
                    #else
                        fvOptions.makeRelative(phiByVoidfraction);
                    #endif
                #else
                    surfaceScalarField phiHbyA
                    (
                        "phiHbyA",
                        (
                            (fvc::interpolate(HbyA) & mesh.Sf() )
                            + phicForces                                   //explicit contribution
                            + rUAfvoidfraction*fvc::ddtCorr(U, phiByVoidfraction)  //correction
                        )
                    );
                    fvOptions.relativeFlux(phiByVoidfraction);
                #endif

                if (modelType=="A")
                    rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                // Update the fixedFluxPressure BCs to ensure flux consistency
                #if defined(version50)
                    if (modelType=="A")
                    {
                        volScalarField rUsed = rUA*voidfraction;
                        constrainPressure(p, U, phiHbyA, rUsed,MRF);
                    }
                    else constrainPressure(p, U, phiHbyA, rUA,MRF);
                #else
                    if (modelType=="A")
                    {
                        setSnGrad<fixedFluxPressureFvPatchScalarField>
                        (
                            #if defined(version40)
                            p.boundaryFieldRef(),
                            #else
                            p.boundaryField(),
                            #endif
                            (
                                phiHbyA.boundaryField()
                                #if defined(version30)
                                - MRF.relative(mesh.Sf().boundaryField() & U.boundaryField())
                                #else
                                - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
                                #endif
                            )/(mesh.magSf().boundaryField()*rUAf.boundaryField()*voidfractionf.boundaryField())
                        );
                    }else
                    {
                        setSnGrad<fixedFluxPressureFvPatchScalarField>
                        (
                            #if defined(version40)
                            p.boundaryFieldRef(),
                            #else
                            p.boundaryField(),
                            #endif
                            (
                                phiHbyA.boundaryField()
                                #if defined(version30)
                                - MRF.relative(mesh.Sf().boundaryField() & U.boundaryField())
                                #else
                                - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
                                #endif
                            )/(mesh.magSf().boundaryField()*rUAf.boundaryField())
                        );
                    }
                #endif

                // Non-orthogonal pressure corrector loop
                #if defined(version30)
                    while (pimple.correctNonOrthogonal())
                #else
                    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                #endif
                {
                    // Pressure corrector
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUAvoidfraction, p) == fvc::div(voidfractionf*phiHbyA) + particleCloud.ddtVoidfraction()
                    );
                    pEqn.setReference(pRefCell, pRefValue);

                    #if defined(version30)
                        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
                        if (pimple.finalNonOrthogonalIter())
                        {
                            phiByVoidfraction = phiHbyA - pEqn.flux()/voidfractionf;
                            phi = voidfractionf*phiByVoidfraction;

                            #include "continuityErrorPhiPU.H"

                            // Explicitly relax pressure for momentum corrector
                            p.relax();

                            U = fvc::reconstruct(phiHbyA)
                              - rUA*fvc::reconstruct(pEqn.flux()/voidfractionf/rUAf);

                            U.correctBoundaryConditions();
                            fvOptions.correct(U);
                        }
                    #else
                        if( corr == nCorr-1 && nonOrth == nNonOrthCorr )
                            #if defined(versionExt32)
                                pEqn.solve(mesh.solutionDict().solver("pFinal"));
                            #else
                                pEqn.solve(mesh.solver("pFinal"));
                            #endif
                        else
                            pEqn.solve();

                        if (nonOrth == nNonOrthCorr)
                        {
                            phiByVoidfraction = phiHbyA - pEqn.flux()/voidfractionf;
                            phi = voidfractionf*phiByVoidfraction;

                            #include "continuityErrorPhiPU.H"

                            // Explicitly relax pressure for momentum corrector
                            p.relax();

                            U = fvc::reconstruct(phiHbyA)
                              - rUA*fvc::reconstruct(pEqn.flux()/voidfractionf/rUAf);

                            U.correctBoundaryConditions();
                            fvOptions.correct(U);
                        }
                    #endif
                } // end non-orthogonal corrector loop

