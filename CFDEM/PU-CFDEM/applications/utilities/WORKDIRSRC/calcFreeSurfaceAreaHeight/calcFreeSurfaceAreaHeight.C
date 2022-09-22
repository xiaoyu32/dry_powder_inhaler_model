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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addRegionOption.H"
    argList::validArgs.append("patchName");
    argList::validArgs.append("direction");
    argList::validArgs.append("charHeight");    
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    const word patchName = args[1];
    const word direction = args[2];
    const word charHeight = args[3];

    char* pEnd;
    scalar height = strtod(charHeight.c_str(), &pEnd);
   
    int idir(0);
    if(direction=="x")
    {
    	idir = 1;
    }
    else if(direction=="y")
    {
    	idir = 2;    	
    }
    else if(direction=="z")
    {
    	idir = 3;    	    
    }
    else
    {
    	FatalError << " Direction is not defined " << abort(FatalError);
    }
    
    Info << "Patch name: " << patchName << endl;
    Info << " " << endl;	
    
    Info << "Calculate surface areas lower the height of " << direction << "=" << height << endl;

    // Create output folder	
    fileName outputRelativePath("FreeSurfaceAreaHeight");
    if( !isDir(mesh.time().path()/outputRelativePath) )
    {
	mkDir(mesh.time().path()/outputRelativePath );
    }

    OFstream* outputFile;
    outputFile =  new OFstream(mesh.time().path()/outputRelativePath/"freeSurfaceAreaHeight.dat");
    *outputFile  << "#Time \t FreeSurfaceArea \t TotalWettedArea" << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();
	
	// Read gas density
	IOobject alphaheader
	(
		"alpha1",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ	
	);

	Info<< " Reading alpha" << endl;
	volScalarField alpha(alphaheader,mesh);	
	
	// Calculate free surface area @ interface
	volScalarField mag_gradAlpha = mag(fvc::grad(alpha));
	
	//scalar totalSurfaceArea = gSum(mag_gradAlpha.internalField()*mesh.V());
        scalar totalSurfaceArea(0);
    	forAll(mesh.cells(),cellI)
	{
		if( mesh.C()[cellI][idir] <= height )
		{
			totalSurfaceArea += mag_gradAlpha.internalField()[cellI]*mesh.V()[cellI];
		}
	}
	// Parallel computation
    	reduce(totalSurfaceArea, sumOp<scalar>());
        
	// Calculate wetted area on the plane
	//const word patchName = "Bottom";
	const label patchI = mesh.boundaryMesh().findPatchID(patchName);	
	
        if (patchI < 0)
        {
        	FatalError
                << "Unable to find patch " << patchName << nl
                << exit(FatalError);
        }
	
	const polyPatch& cPatch = mesh.boundaryMesh()[patchI];

	//scalar totalWettedArea = gSum(alpha.boundaryField()[patchI]*mesh.magSf().boundaryField()[patchI]);
	scalar totalWettedArea(0);
	forAll(cPatch, faceI)
	{
		if( cPatch.faceCentres()[faceI][idir] <= height )
		{
			totalWettedArea += alpha.boundaryField()[patchI][faceI] * mesh.magSf().boundaryField()[patchI][faceI];
		}		
	}	 
	// Parallel computation
    	reduce(totalWettedArea, sumOp<scalar>());
	
    	if(Pstream::master()) //Write only if master
    	{
    		
		Info<< " Writing free surface area into the file " << mesh.time().path()/outputRelativePath/"freeSurfaceAreaHeight.dat" << endl;
		*outputFile	<< alpha.mesh().time().value()  << tab << " " 
                                << totalWettedArea 		<< tab << " " 
				<< totalSurfaceArea 
                		<< endl;

    	}
        
	Info << " " << endl;
    }

    return 0;
}


// ************************************************************************* //
