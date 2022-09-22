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

//===============================
// particle interaction modelling
//===============================

    Info<< "\nReading momentum exchange field Ksl\n" << endl;
    volScalarField Ksl
    (
        IOobject
        (
            "Ksl",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        //dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 1.0)
    );

    Info<< "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "\nCreating dummy density field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh//,
        //dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0), 1.0)
    );
    
    Info<< "Reading particle velocity field Us\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "\nExplicit Source term fsource" << endl;
    volVectorField fsource
    (
        IOobject
        (
            "fsource",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh,
        dimensionedVector("source", dimensionSet(1,-2,-2,0,0),vector(0,0,0)) 
    );
    
//===============================

//#   include "createPhi.H"
#ifndef createPhi_H
#define createPhi_H
Info<< "Reading/calculating face flux field phi\n" << endl;
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

    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");
    bool ignoreFinalPimpleIteration = false;
    if (pimpleDict.found("ignoreFinalPimpleIteration"))
    {
        ignoreFinalPimpleIteration = readBool(pimpleDict.lookup("ignoreFinalPimpleIteration"));
    }
    
    bool momentumCorrectionPeriodicFlow = false;
    scalar rhoParticle(0);
    if (pimpleDict.found("momentumCorrectionPeriodicFlow"))
    {
        momentumCorrectionPeriodicFlow = readBool(pimpleDict.lookup("momentumCorrectionPeriodicFlow"));
	Info << "\n Momentum correction is activated for periodic flow " << endl;
	rhoParticle = readScalar(pimpleDict.lookup("rhoParticle"));
        if (rhoParticle <= 0) FatalError<<" Particle density is not defined in PIMPLE dictionary \n" << abort(FatalError);
    }   

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
 
