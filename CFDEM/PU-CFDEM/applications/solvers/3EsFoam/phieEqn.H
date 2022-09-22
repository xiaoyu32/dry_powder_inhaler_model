// Solve laplacian for potential
solve
(
    fvm::laplacian(phie) + ( alpha1 * charge / volumeP ) / epsilon0
);

// Calculate electrical field
Ef = -fvc::grad(phie);
//Ef.correctBoundaryConditions();

Info<< "Electric field: "
    << "  Min(Ef) = " << min(Ef).value()
    << "  Max(Ef) = " << max(Ef).value()
    << endl; 
