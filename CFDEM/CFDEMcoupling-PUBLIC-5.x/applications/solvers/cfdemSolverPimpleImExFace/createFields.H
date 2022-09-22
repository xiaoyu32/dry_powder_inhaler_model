    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createCFDEMFields.H"

#ifndef createPhi_H
#define createPhi_H
Info<< "Reading/calculating face flux field phi (this is the SUPERFICIAL field, i.e., includes the voidfraction)\n" << endl;
surfaceScalarField phi
 (
     IOobject
     (
        "phi",
        runTime.timeName(),
         mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     linearInterpolate(U*voidfraction) & mesh.Sf()
 );
#endif

 Info<< "Generating interstitial flux field phiByVoidfraction (this is the INTERSTITIAL field!)\n" << endl;
 surfaceScalarField phiByVoidfraction
 (
     IOobject
     (
        "phiByVoidfraction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     linearInterpolate(U) & mesh.Sf()
 );


    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");
    bool ignoreFinalPimpleIteration = false;
    if (pimpleDict.found("ignoreFinalPimpleIteration"))
    {
        ignoreFinalPimpleIteration = readBool(pimpleDict.lookup("ignoreFinalPimpleIteration"));
    }

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );   
