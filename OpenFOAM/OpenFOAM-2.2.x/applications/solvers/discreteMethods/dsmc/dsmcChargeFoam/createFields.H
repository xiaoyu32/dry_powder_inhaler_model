    Info<< "Reading dsmc properties\n" << endl;

    IOdictionary dsmcProperties
    (
        IOobject
        (
            "dsmcProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    dimensionedScalar epsilon0
    (
        dsmcProperties.lookup("epsilon0")
    );

    Info<< "Reading field phie\n" << endl;
    volScalarField phie
    (
	IOobject
	(
            "phie",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
	),
	mesh
    );
    
    Info<< "Creating field Ef\n" << endl;
    volVectorField Ef
    (
        IOobject
        (
            "Ef",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0",dimensionSet(1,1,-3,0,0,-1,0),vector(0,0,0))
    );    

    Info<< "Creating field rhoe\n" << endl;
    volScalarField rhoe
    (
        IOobject
        (
            "rhoe",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ); 
