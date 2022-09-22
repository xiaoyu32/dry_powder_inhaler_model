fvScalarMatrix chargeEqn(charge, charge.dimensions()*dimVol/dimTime);
{   
    electrostaticsKT.compute();   
    
    volScalarField kappaq = electrostaticsKT.kappaq()/rho1;
    if(electrostaticsKT.onTurbulentDiffusivity()) kappaq += electrostaticsKT.kappaqt()/rho1;

    sigmaqEf = electrostaticsKT.sigmaqEf(); 
    sigmaqEf.correctBoundaryConditions();

    chargeEqn =
    (
        (
            fvm::ddt(charge)
          + fvm::div(phi1, charge, "div(phi1,charge)")
          - fvm::Sp(fvc::div(phi1), charge)
        )
      - fvm::laplacian(kappaq, charge)
    );	

    chargeEqn.relax();
    
    //- grad(alpp)/alpp field
    volVectorField normGradAlpp(fvc::grad(alpha1)/fvc::average(alpha1+scalar(1.e-64)));
    
    solve
    (
	chargeEqn == - fvc::div(sigmaqEf/rho1)
      		     - ( normGradAlpp & (sigmaqEf/rho1) )
		     + ( normGradAlpp & (kappaq*fvc::grad(charge)) )	
    );
        
    Info<< "Dispersed phase charge: "
        << "  Min(charge) = " << min(charge).value()
        << "  Max(charge) = " << max(charge).value()
        << endl;    
}
