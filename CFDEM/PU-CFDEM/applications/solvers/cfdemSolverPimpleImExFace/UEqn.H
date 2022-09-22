            surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);

            // Momentum predictor
            #if defined(versionv1606plus) || defined(version40)
                tmp<fvVectorMatrix> tUEqn
            #else
                tmp<fvVectorMatrix> UEqn
            #endif
            (
                fvm::ddt(voidfraction,U)
              + fvm::div(phi, U) //this phi is SUPERFICIAL, i.e., includes the voidfraction!
//              + turbulence->divDevReff(U)
              + particleCloud.divVoidfractionTau(U, voidfraction)
             ==
              - fvm::Sp(Ksl/rho,U)
              + fvOptions(U)
            );

            #if defined(versionv1606plus) || defined(version40)
                fvVectorMatrix& UEqn = tUEqn.ref();
            #endif

            #ifndef version23
                pZones.addResistance(UEqn());
            #endif

            #if defined(versionv1606plus) || defined(version40)
                UEqn.relax();
                fvOptions.constrain(UEqn);
                volScalarField rUA = 1.0/UEqn.A();
            #else
                UEqn().relax();
                fvOptions.constrain(UEqn());
                volScalarField rUA = 1.0/UEqn().A();
            #endif

            surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

            surfaceScalarField phicForces
            (
                fvc::interpolate(rUA*(Ksl*Us - f)/rho) & mesh.Sf()
            );

            #if defined(version30)
                if (pimple.momentumPredictor())
            #else
                if (momentumPredictor)
            #endif
            {
                if (modelType=="B" || modelType=="Bfull")
                    #if defined(versionv1606plus) || defined(version40)
                        solve(
                            UEqn == fvc::reconstruct(phicForces/rUAf
                                      - fvc::snGrad(p)*mesh.magSf())
                        );
                    #else
                        solve(
                            UEqn() == fvc::reconstruct(phicForces/rUAf
                                        - fvc::snGrad(p)*mesh.magSf())
                        );
                    #endif
                else
                    #if defined(versionv1606plus) || defined(version40)
                        solve(
                            UEqn == fvc::reconstruct(phicForces/rUAf
                                      - fvc::snGrad(p)*voidfractionf*mesh.magSf())
                        );
                    #else
                        solve(
                            UEqn() == fvc::reconstruct(phicForces/rUAf
                                        - fvc::snGrad(p)*voidfractionf*mesh.magSf())
                        );
                    #endif

                fvOptions.correct(U);
            }
