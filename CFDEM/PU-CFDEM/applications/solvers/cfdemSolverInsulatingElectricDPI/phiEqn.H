Info<< "Entering phiEqn" << endl;

//Solve Poisson equation for electric potential

for (int phieCorr = 0; phieCorr<=nPhieCorr; phieCorr++)
{
    // -- solve electric potential for the interior domain --
    fvScalarMatrix phieEqn
    (
	fvm::laplacian(mixturePermittivity,phie) + rhoe/epsilon0
    );

    //label phieRefCell = 0;
    //scalar phieRefValue = 0.0;
    //phieEqn.setReference( phieRefCell, phieRefValue );
    
    phieEqn.solve();
    
    // -- solve electric potential for the exterior domain --
    fvScalarMatrix envPhieEqn( 
				    fvm::laplacian(envPhie) 
			    );
    envPhieEqn.solve();    
    
}

phie.correctBoundaryConditions();

//Solve electrical field
Ef = -fvc::grad(phie);

//Update electric field gradient
gradEf = fvc::grad(Ef);

Info<< "Leaving phiEqn" << endl;


