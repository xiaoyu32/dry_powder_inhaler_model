
vector fluidMomentum = fvc::domainIntegrate(rho*voidfraction*U).value();
vector solidMomentum = fvc::domainIntegrate(rhoParticle*(1-voidfraction)*Us).value();
scalar domainRhomix  = fvc::domainIntegrate(rhoParticle*(1-voidfraction)).value()+fvc::domainIntegrate(rho*voidfraction).value();

volVectorField momentumCorrection
(
    IOobject
    (
        "momCorr",
        particleCloud.mesh().time().timeName(),
        particleCloud.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    particleCloud.mesh(),
    dimensionedVector( "zero", dimensionSet(0,1,-2,0,0),
    ( fluidMomentum + solidMomentum ) / particleCloud.mesh().time().deltaT().value() / domainRhomix )
);

Info << " " << endl;
Info << "Domain-averaged fluid momentum = " << fluidMomentum << endl;
Info << "Total particle momentum = " << solidMomentum << endl;
Info << "Total momentum = " << fluidMomentum + solidMomentum << endl;
