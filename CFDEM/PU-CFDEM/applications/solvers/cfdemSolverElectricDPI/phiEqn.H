Info<< "Entering phiEqn" << endl;
{
//Solve Poisson equation for electric potential

fvScalarMatrix phieEqn
(
    fvm::laplacian(mixturePermittivity,phie) + rhoe/epsilon0
);

label phieRefCell = 0;
scalar phieRefValue = 0.0;

phieEqn.setReference( phieRefCell, phieRefValue );
phieEqn.solve();

//Solve electrical field
Ef = -fvc::grad(phie);

//Update electric field gradient
gradEf = fvc::grad(Ef);

}
Info<< "Leaving phiEqn" << endl;


