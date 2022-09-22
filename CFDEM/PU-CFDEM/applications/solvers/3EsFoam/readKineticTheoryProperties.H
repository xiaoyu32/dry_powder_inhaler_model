    IOdictionary kineticTheoryProperties
    (
        IOobject
        (
            "kineticTheoryProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar preAlphaExp(readScalar(kineticTheoryProperties.lookup("preAlphaExp")));

    dimensionedScalar alphaMax(kineticTheoryProperties.lookup("alphaMax"));

    scalar expMax(readScalar(kineticTheoryProperties.lookup("expMax")));

    dimensionedScalar g0(kineticTheoryProperties.lookup("g0"));

    Switch packingLimiter(kineticTheoryProperties.lookup("packingLimiter"));

    dimensionedScalar alphaRes(kineticTheoryProperties.lookup("alphaResidual"));
