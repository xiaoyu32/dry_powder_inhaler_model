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

#ifndef domainFilter_H
#define domainFilter_H

#include "fluidParticleCoarsening.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
	void createStencils
	(
		const fvMesh& mesh,
		const int& filterwidth_,
		const labelListList& mA_,
		const int& max_Ax_,
		const int& max_Ay_,
		const int& max_Az_,
		labelListList& stencillist_
	);

	void constructfilter
	(
	        const argList& args,
	        const Time& runTime,
	        const fvMesh& mesh,
	        const int& minfilterwidth,
	        const int& maxfilterwidth,
	        const int& filterincrement,
	        labelListListList& stencillistfilter,
	        const fileName& outputRelativePath_
	);
}	

#endif