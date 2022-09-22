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

#include "calcCylindEf.H"
#include "timeSelector.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include <fstream>
#include <iostream>
#include "math.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
#   include "addRegionOption.H"
    Foam::argList::addBoolOption
    (
        "noWrite",
        "suppress writing results"
    );
    #include "addDictOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"   
    #include "createFields.H"
    
    
    
    forAll(timeDirs, timeI)
    {

	runTime.setTime(timeDirs[timeI], timeI);
	if ( runTime.timeName() != "0" )
	{
        					
		Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;

        	mesh.readUpdate();
		
		// Read electric field
		
        	IOobject Efheader
        	(
            		"Ef",
            		runTime.timeName(),
            		mesh,
            		IOobject::MUST_READ
        	);
		
		volVectorField Ef_f(Efheader,mesh);
		
		label nCells = returnReduce(mesh.cells().size(), sumOp<label>());
		
		
		std::ofstream writer;
	
		writer.open("EfData1.dat");
		
		writer<<nCells<<nl;
		
		
		forAll ( Ef_f, iCell)
		{
			Foam::vector cellEf = Ef_f[iCell];
			Foam::vector cellPos = Ef_f.mesh().C()[iCell];
			
			writer<< cellEf.x() <<" ";
			writer<< cellEf.y() <<" ";
			writer<< cellEf.z() <<" "; 
			writer<< Foam::sqrt(cellEf.x()*cellEf.x()+cellEf.y()*cellEf.y()) <<" "; // r
			writer<< Foam::atan2(cellEf.y(),cellEf.x()) <<" "; // theta
			
			writer<< cellPos.x() <<" ";
			writer<< cellPos.y() <<" ";
			writer<< cellPos.z() <<" ";  // z	
			writer<< Foam::sqrt(cellPos.x()*cellPos.x()+cellPos.y()*cellPos.y()) <<" "; // r
			writer<< Foam::atan2(cellPos.y(),cellEf.x()) <<nl; // theta	 
		}
		
		writer.close();       	
				
	}
	
		
    }
    

	
    return 0;
}


// ************************************************************************* //
