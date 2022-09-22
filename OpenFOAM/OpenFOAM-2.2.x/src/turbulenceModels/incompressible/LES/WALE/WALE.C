#include "WALE.H"
#include "addToRunTimeSelectionTable.H"
//#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WALE, 0);
addToRunTimeSelectionTable(LESModel, WALE, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WALE::WALE
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenEddyVisc(U, phi, transport),

    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    ),
   cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    ),
  ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce",
            coeffDict_,
            1.048
        )
    ),
 nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )


{

    printCoeffs();
}


tmp<volSymmTensorField>WALE::Sd
(
    const volTensorField& gradU
) const
{
    return dev(symm(gradU & gradU));
}


tmp<volScalarField> WALE::k
(
    const volTensorField& gradU
) const
{
    volScalarField magSqrSd(magSqr(Sd(gradU)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                runTime_.timeName(),
                mesh_
            ),
            sqr(sqr(cw_)*delta()/ck_)*
            (
                pow3(magSqrSd)
               /(
                   sqr
                   (
                       pow(magSqr(symm(gradU)), 5.0/2.0)
                     + pow(magSqrSd, 5.0/4.0)
                   )
                 + dimensionedScalar
                   (
                       "SMALL",
                       dimensionSet(0, 0, -10, 0, 0),
                       SMALL
                   )
               )
            )
        )
    );
}


void WALE::correctNut()
{
    nut_ = ck_*delta()*sqrt(k(fvc::grad(U_)));
    nut_.correctBoundaryConditions();
    //fv::option::New(this->mesh_).correct(this->nut_);
   //  fvOptions.correct(nut_);
    correctNut();

}


tmp<volScalarField> WALE::epsilon() const
{
    volScalarField k(this->k(fvc::grad(U_)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ce_*k*sqrt(k)/delta()
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WALE::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);

}


bool WALE::read()
{
    if (GenEddyVisc::read())
    {
        ck_.readIfPresent(coeffDict());
        cw_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
