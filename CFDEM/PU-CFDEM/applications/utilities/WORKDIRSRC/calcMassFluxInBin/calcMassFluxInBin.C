/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    postCalc

Description
    Generic wrapper for calculating a quantity at each time

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"

#include "timeSelector.H"
#include "fvCFD.H"

#include "IOmanip.H"
#include "OFstream.H"

//#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addRegionOption.H"
    argList::validArgs.append("direction");
    argList::validArgs.append("charHeight");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    const word direction = args[1];
    const word charHeight = args[2];

    int idir(0);
    if(direction=="x")
    {
        idir = 0;
    }
    else if(direction=="y")
    {
        idir = 1;
    }
    else if(direction=="z")
    {
        idir = 2;
    }
    else
    {
        FatalError << " Direction is not defined " << abort(FatalError);
    }

    char* pEnd;
    scalar height = strtod(charHeight.c_str(), &pEnd);

    // Create output file 

    OFstream* outputFile;
    outputFile =  new OFstream
    (
   	runTime.processorCase()
 	? runTime.path()/".."/"solidInBin.dat" // DEMinputFilename
 	: runTime.path()/"solidInBin.dat"      // DEMinputFilename
    );
    *outputFile  << "#Time \t averagedSolidVolFrac \t volumetricSolidFlux" << endl;

    scalar weightedAlphaPre(0);
    scalar weightedAlphaIns(0);	
    scalar diffWeightedAlpha(0);

    /*
    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<phaseModel> phasea = phaseModel::New
    (
        mesh,
        transportProperties,
        "a"
    );
      
    const dimensionedScalar& rhoa = phasea->rho();
    */	
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();
	
	// Read gas density
	IOobject alphaheader
	(
		"alpha",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ	
	);

	//Info<< " Reading alpha" << endl;
	volScalarField alpha(alphaheader,mesh);	
	
	scalar weightedAlpha(0);
	scalar totalVol(0);
	forAll(mesh.cells(),cellI)
	{
		if( mesh.C()[cellI][idir] > height)
		{
			totalVol += mesh.V()[cellI];
			weightedAlpha += alpha[cellI] * mesh.V()[cellI];
		}
	}

	weightedAlphaIns = weightedAlpha; 
	scalar deltaT(0);

	if(timeI==0)
	{
		deltaT = 1;
		diffWeightedAlpha = weightedAlphaIns;
	}else
	{
		deltaT = timeDirs[timeI].value() - timeDirs[timeI-1].value();
		diffWeightedAlpha = (weightedAlphaIns - weightedAlphaPre)/deltaT;
	}

        if(Pstream::parRun())
	{
		reduce(diffWeightedAlpha, sumOp<scalar>());
		reduce(totalVol, sumOp<scalar>());
	}          

	if(Pstream::master())
	{
        	Info << " Averaged solid vol. frac. = " << weightedAlphaIns/totalVol << endl;
        	Info << " Volumetric solid flux = " << diffWeightedAlpha << endl;

		if(timeDirs[timeI].value() != 0)
		{
			*outputFile	<< alpha.mesh().time().value()  << tab << " " 
                		        << weightedAlphaIns/totalVol	<< tab << " " 
					<< diffWeightedAlpha		<< tab << " " 
                			<< endl;
        	}
	}

	weightedAlphaPre = weightedAlphaIns; 

	Info << " " << endl;
    }

    return 0;
}


// ************************************************************************* //
