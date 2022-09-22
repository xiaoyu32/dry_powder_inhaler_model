{
    
    //fluxes
    surfaceScalarField uflux = fvc::interpolate( voidfraction ) * phi;
    carrierMassFlux = linearInterpolate(voidfraction) * phi;
    
    //advance passive scalar fields
    for( label ii = 0; ii < passiveScalar.size(); ++ii )
    {

	fvScalarMatrix passiveScalarEqn
	(
            fvm::ddt( voidfraction, passiveScalar[ii]() ) 
	  + fvm::div( uflux, passiveScalar[ii]() )
            ==
      	    fvm::laplacian( kPassiveScalar * voidfraction, passiveScalar[ii]() )
	  + fvm::Sp( impPassiveScalarSource[ii](), passiveScalar[ii]() ) 
	  + expPassiveScalarSource[ii]()
	);
            
	passiveScalarEqn.relax();
	passiveScalarEqn.solve();

	//update passive scalar flux
    	passiveScalarFlux[ii]() = carrierMassFlux * linearInterpolate( passiveScalar[ii]() );
	
    }
    
}

