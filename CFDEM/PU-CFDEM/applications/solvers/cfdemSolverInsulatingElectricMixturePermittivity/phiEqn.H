Info<< "Entering phiEqn" << endl;

// XL: calculate mixture permittivity from Bruggeman equation

forAll( mixturePermittivity, iCell )
{
   
   Foam::dimensionedScalar cellPermittivity = epsilonGas * voidfraction[iCell] + epsilonSolid * ( 1.0 - voidfraction[iCell] );
   
   
   scalar cellPermittivity_ = cellPermittivity.value();
   const scalar epsilonSolid_ = epsilonSolid.value();
   const scalar epsilonGas_ = epsilonGas.value();
   const scalar alphacValue = voidfraction[iCell];
      
   for( int ii = 0; ii < 3; ++ii )
   {
       const scalar ff = -cellPermittivity_ +  
       			(  
       			     pow
			     ( 
       		       		( epsilonSolid_ - cellPermittivity_ )/
				( epsilonSolid_ - epsilonGas_ ), 3
			     ) 
			     * epsilonGas_/pow( alphacValue, 3 ) 
			);

       const scalar gradff = -1.0 + 
       			     (   -3.0 *
       				  pow
				  ( 
       		       		     ( epsilonSolid_ - cellPermittivity_ ), 2
				  )	
				  * 
				  pow
				  ( 
       		       		     ( epsilonSolid_ - epsilonGas_ ), 3
				  )
				  * epsilonGas_/pow( alphacValue, 3 )
			     );			
			
       cellPermittivity_ = cellPermittivity_ - ff/gradff;
   }
   
   mixturePermittivity[iCell] = cellPermittivity_;
   
}

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


