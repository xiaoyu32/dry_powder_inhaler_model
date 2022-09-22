//- Update charge density
rhoe.internalField() = dsmc.rhoe();

//-Solve laplacian for potential
solve
(
    fvm::laplacian(phie) + rhoe / epsilon0
);

// Calculate electrical field
Ef = -fvc::grad(phie);

Info<< "Electric field: "
    << "  Min(Ef) = " << min(Ef).value()
    << "  Max(Ef) = " << max(Ef).value()
    << endl; 
