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

#ifndef filteringEulerianField_H
#define filteringEulerianField_H

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

#include "./ann/include/ANN/ANN.h"

#include "CPCCellToCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{	
	void constructfilter
	(
		const argList& args, 
		const Time& runTime, 
		const fvMesh& mesh,
		const int& minfilterwidth,
		const int& maxfilterwidth,
		const int& filterincrement,
		labelListListList& stencillistfilter
	);
	
	void filteringEulerianVariables
	(
		const argList& args, 
		const Time& runTime, 
		const fvMesh& mesh,
		const labelListListList& stencillistfilter,
		const int& filterwidth,
		const volScalarField& voidfraction_,
		const volVectorField& U_,
		volScalarField& alpf_,
		volVectorField& tildeUf_ 
	);
	
	void DragForce
	(
		cfdemCloud& sm,
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volScalarField& rho_,
		const bool& verbose_,
		vectorField& DragForce_
	);		


	void filteredDragCoefficient
	(
		cfdemCloud& sm,
		const bool& verbose_,
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volVectorField& tildeUf_,
		const volScalarField& barPg_,
		const volScalarField& p_,
		const volScalarField& rho_,
		const vectorField& DragForce_,
		vectorField& ResolvedDragForce_,
		vectorField& gii_
	);
	
 	void ParticleCoarsening
	(
		const fvMesh& mesh,
		cfdemCloud& sm,
		const int& nparticle_,
		const bool verbose_,
		const volVectorField& Uf_,
		const volVectorField& tildeUf_,
		const volScalarField& p_,
		const volScalarField& barPg_,
		const vectorField& DragForce_,
		const vectorField& ResolvedDragForce_,
		vectorField& parcelgii_			
	);
	
	void applyBins
	(
		cfdemCloud& sm,
		const int& filterwidth,
		const int& maxfilterwidth,
		const int& nBin_,
		const scalar& maxalpp_,
		const volScalarField& baralpf_,
		const volVectorField tildeUf_,
		vectorField& gii_,
		Field< Field < Field <scalar> > >& giiCondalppFilter_,
		Field< Field < Field <scalar> > >& numbergiiCondalppFilter_
	);

	void writeBins
	(
	    	const fvMesh& mesh,
		const int& nBin_,
		const scalar& maxalpp_,
		const fileName outputRelativePath_,		
		const int& minfilterwidth,
		const int& maxfilterwidth,
		const int& filterincrement,
		const int& nparticle_,
		Field< Field < Field <scalar> > >& giiCondalppFilter_,
		Field< Field < Field <scalar> > >& numbergiiCondalppFilter_			
	);				
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
