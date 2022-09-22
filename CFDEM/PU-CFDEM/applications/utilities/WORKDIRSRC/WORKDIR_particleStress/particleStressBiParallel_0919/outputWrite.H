/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Function prototype for all simple post-processing functions
    e.g. filterinEulerianFieldDivPhi, filterinEulerianFieldMagU etc.

\*---------------------------------------------------------------------------*/

#ifdef outputWrite_H
#define outputWrite_H

#include "particleStressBiParallel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{	
	// Initiate output file
	void initOutput
	(	
		const int&          nParticleClass_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		OFstream*&          outputCollFile_,
		const bool&             bidisperse_,
		const bool&          calcCollision_		);
	
	// Write into files	
	void writeOutput
	(	
		const int&          nParticleClass_,
		const Time& 		   runTime_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		OFstream*&          outputCollFile_,
		const bool&             bidisperse_,
		const bool&          calcCollision_,
		SymmTensor<double>*       sigmaKin_,					
		SymmTensor<double>*	 sigmaColl_,
		const scalar&		 domainVol_,
		const scalar&	   alppInSubVolume_		);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
