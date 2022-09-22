    // Mean velocities in each cell
    volVectorField meanVelCell
    (
	 IOobject
	 (
            "meanUp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	  ),
	mesh,
	dimensionedVector( "zero", dimensionSet(0,1,-1,0,0), vector(0,0,0) )
    );

    volScalarField mappedVoidfraction
    (
        IOobject
        (
            "mappedVoidfraction",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0),scalar(1.)) 
    );

    volVectorField mappedUs
    (
        IOobject
        (
            "mappedUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0),vector(0,0,0)) 
    );

    volScalarField weightField
    (
        IOobject
        (
            "weightField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0),scalar(0)) 
    );
