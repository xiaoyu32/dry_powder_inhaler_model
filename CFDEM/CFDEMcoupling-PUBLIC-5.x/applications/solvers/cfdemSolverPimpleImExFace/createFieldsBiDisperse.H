    Info<< "\nGeneration Additional Fields relevant for Bi-disperse Calculations \n" << endl;
    volScalarField phiP1
    (
        IOobject
        (
            "phiP1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        voidfraction
    );

    volScalarField phiP2
    (
        IOobject
        (
            "phiP2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        voidfraction
    );

    volScalarField dSauter
    (
        IOobject
        (
            "dSauter",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,1,0,0,0), 0.0)
    );

    volVectorField uP1
    (
        IOobject
        (
            "uP1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        U
    );

    volVectorField uP2
    (
        IOobject
        (
            "uP2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        U
    );
