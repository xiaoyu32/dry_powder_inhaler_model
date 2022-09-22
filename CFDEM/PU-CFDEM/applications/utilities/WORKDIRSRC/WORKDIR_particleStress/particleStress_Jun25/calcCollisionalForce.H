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

#ifndef calcCollisionalForce_H
#define calcCollisionalForce_H

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

#include "./ann/include/ANN/ANN.h"

#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{	

	void calcLiqIn(double*& liq, int& i, int& j, double**& liquidVol);
	
	void calcLiqBack(double *&liq, int& i, int& j, double**& liquidVol);

	void calcForce(cfdemCloud& particles, 	int collisionModel,
						int j, 		// Particle	
						int i,  	// neighboring particles
				        	double **& fcoll, 
						double **& ftan,
						double& k_n, 
						double& k_t, 
						double& gamma_n, 
						double& gamma_t,
						double& youngsModulus,		  
						double& poissonsRatio,
						double& coefficientRestitution,		
						double& coefficientFriction,					 
						double***& delta_t, 
						double& mu_f, 
						const scalar& dt, 
						bool& tangential_history, 
						double *&liq, 
						bool& liquid_transfer, 
						double **& liquidVol, 
						double& surf_tension, 
						double& fluid_visc, 
						double **& fcap, 
						double **& fvisc,
						bool**& first_touch,
						bool& cohesion,
						double& minimumDistanceVdW,
						double& cohEnergyDens,
						double **& fcoh,
						const scalar& rhop,
						double& e_n,
						double& e_t,
						symmTensor& simga_coll	   	    );

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
