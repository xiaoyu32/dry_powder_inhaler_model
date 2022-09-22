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

#ifndef createParcels_H
#define createParcels_H

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"

#include "cfdemCloud.H"
#include "dataExchangeModel.H"
#include "voidFractionModel.H"
#include "locateModel.H"
#include "averagingModel.H"
#include "momCoupleModel.H"
#include "forceModel.H"
#include "IOModel.H"
#include "interpolationCellPoint.H"

#include "timeSelector.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "../ann/include/ANN/ANN.h"

#include "IOmanip.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{	

 	void createParcels
	(
		const fvMesh& mesh,
		ANNkd_tree*& kdTree,			
		cfdemCloud& sm,
		const int& parcelSize_,		
		int**& parcelCloud_,
		double**& parcelPositions_,
		double**& parcelVelocities_,
		int*& parcelNparts_,				
		double** & parcelKinStress_,
		scalar& aveSubQparcel2_,
		vector& meanParcelVel_,		
		const bool verbose_ 			
	);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
