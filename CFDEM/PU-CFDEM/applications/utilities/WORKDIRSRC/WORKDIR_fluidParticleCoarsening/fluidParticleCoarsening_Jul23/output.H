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

#ifndef output_H
#define output_H

#include "fluidParticleCoarsening.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
	void writeBins
	(
		const Time& runTime,
	    	const fvMesh& mesh,
		const int& nBin_,
		const scalar& minalpp_,
		const scalar& maxalpp_,
		const scalar& minVelSlip_,
		const scalar& maxVelSlip_,
		const scalar& minTaupVelSlip_,
		const scalar& maxTaupVelSlip_,
		const fileName outputRelativePath_,		
		const int& minfilterwidth,
		const int& maxfilterwidth,
		const int& filterincrement,
		const int& nparticle_,
		Field< Field < Field <scalar> > >& 	giiCondalppFilter_,
		Field< Field < Field <scalar> > >& 	velSlipCondalppFilter_,
		Field < Field <scalar> > & 		taupVelSlipCondalppFilter_,
		Field < Field <scalar> > & 		numberCondalppFilter_,	
		Field< Field < Field <scalar> > >&	VelSlipJointgiiCondalppFilter_,
		Field< Field < Field <scalar> > >&	TaupVelSlipJointgiiCondalppFilter_,	
		Field <Field < Field <scalar> > >& 	numberVelSlipJointgiiCondalppFilter_,
		Field <Field < Field <scalar> > >& 	numberTaupVelSlipJointgiiCondalppFilter_,				
		bool verbose_,
		bool EulerEulerFiltering_,
		bool EulerianVslipBin_ 
	);			
}	

#endif