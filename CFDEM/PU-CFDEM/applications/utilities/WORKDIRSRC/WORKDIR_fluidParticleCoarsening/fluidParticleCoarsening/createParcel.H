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

#ifndef createParcel_H
#define createParcel_H

#include "fluidParticleCoarsening.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
 	void createParcels
	(
		const fvMesh& mesh,			
		cfdemCloud& sm,
		cfdemCloud& parcelsm,
		const int& nparticle_,
		const bool verbose_,
		const fileName outputRelativePath_, // Write parcel diameter in a file in the directory of ./PostProcessing..	
		labelListList& parcelList,
		const bool& weighting_,
		const scalar& minX_,
		const scalar& maxX_,
		const scalar& minY_,
		const scalar& maxY_,
		const scalar& minZ_,
		const scalar& maxZ_,
	    //ANNpointArray& dataPts,
	    //ANNpoint& queryPt,
	    //ANNkd_tree*& kdTree,
		const labelList& partIDInSubVolume_,
		labelList& parcelIDInSubVolume_		 			
	);

}	

#endif
