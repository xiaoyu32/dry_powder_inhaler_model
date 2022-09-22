      particleCloud.clockM().start(2,"Coupling");
      bool hasEvolved =
            particleCloud.evolve(voidfraction,Us,U); //returns already smoothened fields!

      #if defined(biDisperseModels)
      if(doBiDisperseUpdates)  particleCloud.updateEulerArrays(voidfraction, phiP1, phiP2, dSauter, uP1,  uP2);
      #endif

      //Repartition and Smoothen Coupling Forces if necessary
      if(hasEvolved)
      {
        //Repartition the implicit and explicit forces acting on the fluidVel
        particleCloud.forceM(0).repartitionImExForces();

        //Smoothen implicit momCoupling force
        particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());

        //Smoothen explicit momCoupling force
        particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).expParticleForces());

        //Handle regions void of particles. Must be done after smoothing
        particleCloud.forceM(0).treatVoidCells();

      }

      //Compute explicit momCoupling if there is any - consider time step fraction, so do every fluid time step!
      label idExpCouple(particleCloud.registryM().getProperty("explicitCouple_index"));
      if(idExpCouple>=0)
      {
            #if defined(version40)
            f
               =  particleCloud.momCoupleM(particleCloud.registryM().getProperty("explicitCouple_index")).expMomSource();
            #else
            f.internalField()
               =  particleCloud.momCoupleM(particleCloud.registryM().getProperty("explicitCouple_index")).expMomSource();
            #endif
            f.relax();
            f.correctBoundaryConditions();
      }

      //Compute coupling strength - consider time step fraction, so do every fluid time step!
      #if defined(version40)
          Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
      #else
          Ksl.internalField() = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
      #endif
      Ksl.relax();
      Ksl.correctBoundaryConditions();

      //Force Checks - will set fImpTotal and fExpTotal
      #include "forceCheckImEx.H"
      #include "couplingForceError.H"

      particleCloud.clockM().stop("Coupling");
      particleCloud.clockM().start(10,"Flow");

