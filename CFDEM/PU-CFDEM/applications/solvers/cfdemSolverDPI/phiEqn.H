Info<< "Entering phiEqn" << endl;
{
//Solve Poisson equation for electric potential

fvScalarMatrix phieEqn
(
    fvm::laplacian(phie) == rhoe/epsilon0
);

label phieRefCell = 0;
scalar phieRefValue = 0.0;

phieEqn.setReference( phieRefCell, phieRefValue );
phieEqn.solve();

//Solve electrical field
Ef = fvc::grad(phie);
divE = epsilon0 * fvc::div(Ef);//epsilon0 * fvc::laplacian(phie);//


//Ef = fvc::reconstruct(fvc::snGrad(phie) * mesh.magSf());

}
Info<< "Leaving phiEqn" << endl;


