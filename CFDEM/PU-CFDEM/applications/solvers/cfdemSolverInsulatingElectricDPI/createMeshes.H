// The primary mesh for the flow field
Foam::Info
    << "Create mesh for the flow region for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);

// The surrounding region for calculating the electric field

Foam::IOdictionary regionDict = Foam::IOdictionary
(
    IOobject
    (
        "regionProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT
    )
);


const word envRegionName = regionDict.lookupOrDefault("envRegionName", word("env"));

Foam::Info
	<< "Create mesh for region '" << envRegionName << "'." << Foam::endl;

Foam::fvMesh envMesh
(
	Foam::IOobject
	(
		envRegionName,
		runTime.timeName(),
		runTime,
		Foam::IOobject::MUST_READ
	)
);




