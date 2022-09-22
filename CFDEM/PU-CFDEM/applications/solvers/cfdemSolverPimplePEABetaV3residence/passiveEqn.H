{
    
    /*fvScalarMatrix passiveScalarEqn
    (
        fvm::ddt(voidfraction, passiveScalar) 
      + fvc::div
        (
                voidfraction * U * passiveScalar
        )
	==
        kPassiveScalar * fvm::laplacian(passiveScalar)
    );*/
    
    //surfaceScalarField uflux = fvc::interpolate( voidfraction * U ) & mesh.Sf();
    surfaceScalarField uflux = fvc::interpolate( voidfraction ) * phi;
    //surfaceScalarField uflux = fvc::interpolate( U ) & mesh.Sf();
    
    /*fvScalarMatrix passiveScalarEqn
    (
        fvm::ddt( passiveScalar ) 
      + fvm::div( uflux, passiveScalar )
        ==
      	fvm::laplacian( kPassiveScalar, passiveScalar )
    );*/
    
    fvScalarMatrix passiveScalarEqn
    (
        fvm::ddt( voidfraction, passiveScalar ) 
      + fvm::div( uflux, passiveScalar )
        ==
      	fvm::laplacian( kPassiveScalar * voidfraction, passiveScalar )
    );
        
    passiveScalarEqn.relax();
    passiveScalarEqn.solve();
    
    //update fluxes
    carrierMassFlux = linearInterpolate(voidfraction) * phi;
    passiveScalarFlux = carrierMassFlux * linearInterpolate(passiveScalar);
    
}

