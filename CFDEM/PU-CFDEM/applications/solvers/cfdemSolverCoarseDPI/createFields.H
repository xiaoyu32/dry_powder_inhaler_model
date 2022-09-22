    Info<< "Reading physicalProperties\n" << endl;
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
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

    Info<< "\nResidual source term fresidual" << endl;
    volVectorField fresidual
    (
        IOobject
        (
            "fresidual",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedVector("residual", dimensionSet(1,-2,-2,0,0),vector(0,0,0)) 
    );
   
   
    volScalarField apiWallDensityVol
    (
         IOobject
        (
            "apiWallDensityVol",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedScalar( "zero", dimensionSet(0,-2,0,0,0),scalar(0) )   
    );
   
        
//===============================

//#   include "createPhi.H"
//#ifndef createPhi_H
//#define createPhi_H

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
     linearInterpolate(U) & mesh.Sf()
 );
//#endif

    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");
    
    scalar rhoParticle(0);	
    bool fullyPeriodicFLow(false);

    bool mixtureVelocityCorrectionPeriodicFlow(false);
    if (pimpleDict.found("mixtureVelocityCorrectionPeriodicFlow"))
    {
        mixtureVelocityCorrectionPeriodicFlow = true;
        Info << "Mixture velocity correction is activated for periodic flow " << endl;
    } 
    
    // Boundary patches, find coupled patches
    label nCoupledPatches(0);
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].type() == "cyclic" )
        {
        	nCoupledPatches++;
	}
    }

    //if ( nCoupledPatches == 6 ) 
   if ( nCoupledPatches > 2 )  
   {
	Info << "Number of coupled patches = " << nCoupledPatches << endl;
        fullyPeriodicFLow = true; 
	Info << "Fully periodic flow " << endl;
	rhoParticle = readScalar(pimpleDict.lookup("rhoParticle"));
    }	

    bool correctCouplingError(false);
    if (pimpleDict.found("correctCouplingError"))
    {
        correctCouplingError = true;
	Info << "Momentum coupling error correction is active " << endl;
    }
   
    //if ( nCoupledPatches < 6 )
    if ( nCoupledPatches < 2 )
    {
	mixtureVelocityCorrectionPeriodicFlow = false;
    }	
     
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
    
    //- 
    bool calcStress(false);
    if (pimpleDict.found("calculateParticleStresses"))
    {
        calcStress = true;
	Info << "Solver calls DEM solver to compute particle stresses " << endl;
    }    
 
