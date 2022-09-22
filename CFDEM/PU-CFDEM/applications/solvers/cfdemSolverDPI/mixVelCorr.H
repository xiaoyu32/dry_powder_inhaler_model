vector fluidMomentum = fvc::domainIntegrate(rho*voidfraction*U).value();
vector solidMomentum = fvc::domainIntegrate(rhoParticle*(1-voidfraction)*Us).value();

volVectorField velocityCorrection
(   
	IOobject
	(
	    "velCorr",
	    particleCloud.mesh().time().timeName(),
	    particleCloud.mesh(),
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	particleCloud.mesh(),
	dimensionedVector( "zero", dimensionSet(0,1,-1,0,0), 
	( fluidMomentum + solidMomentum ) / ( fvc::domainIntegrate(rhoParticle*(1-voidfraction)).value()+fvc::domainIntegrate(rho*voidfraction).value() ) )
);

Info << "Mixture velocity correction " << endl;
Info << "Velocity correction at cell[0] = " << velocityCorrection[0] << endl;

// Total momentum correction for fluid and solid
U -= velocityCorrection;
Us-= velocityCorrection; 
 
 
