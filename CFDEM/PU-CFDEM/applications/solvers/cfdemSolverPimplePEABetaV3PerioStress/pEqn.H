{
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           (fvc::interpolate(HbyA) & mesh.Sf())
//         + voidfractionf*rAUf*fvc::ddtCorr(U, phic)
         + fvc::ddtPhiCorr(rAU*voidfraction, U, phi)
         + phiForces
        )
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
//    setSnGrad<fixedFluxPressureFvPatchScalarField>
//    (
//        p.boundaryField(),
//        (
//            phiHbyA.boundaryField()
//            - (mesh.Sf().boundaryField() & U.boundaryField())
//        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
//    );

    // Fix for gravity on outlet boundary.
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiForces.boundaryField()[patchi] = 0.0;
        }
    }

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(voidfractionf*rAUf, p)
         ==
            fvc::ddt(voidfraction) + fvc::div(voidfractionf*phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux()/voidfractionf;

            p.relax();

            U = HbyA
              + rAU*fvc::reconstruct((phiForces - pEqn.flux()/voidfractionf)/rAUf);
            U.correctBoundaryConditions();
	    fvOptions.correct(U);
        }
    }
}

#include "continuityErrs.H"
