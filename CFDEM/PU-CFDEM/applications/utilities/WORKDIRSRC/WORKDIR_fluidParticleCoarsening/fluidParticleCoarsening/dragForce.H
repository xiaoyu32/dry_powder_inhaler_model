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

#ifndef dragForce_H
#define dragForce_H

#include "fluidParticleCoarsening.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
	void WenYuDragForce
	(		
		const vector& Ur_,
		const scalar& ds_,
		const scalar& rho_,
		const scalar& nuf_,
		const scalar& voidfraction_,
		vector& WenYuDrag_
	);
		
	void CalculateDragForce
	(
		cfdemCloud& sm,
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volScalarField& rho_,
		const bool& verbose_,
		vectorField& DragForce_,
		const labelList& partIDInSubVolume_
	);
		
	void EulerianParticleVelocityForce
	(
		cfdemCloud& sm,			
		const fvMesh& mesh,
		volVectorField&	Up_,
		const vectorField& DragForce_,
		volVectorField& MappedDragForce_,
		const labelList& partIDInSubVolume_			
	);

	void CalculateEulerianDragForce
	(		
		cfdemCloud& sm,
		const fvMesh& mesh,
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volVectorField& Up_,
		const volScalarField& rho_,
		volVectorField& EulerianDragForce_	
	);	
	
	void filteredEulerEulerDragCoefficient
	(
		cfdemCloud& sm,
		const bool& verbose_,
		const fvMesh& mesh,
		const labelListListList& stencillistfilter,
		const int& filterwidth,		
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volVectorField& Up_,		
		const volScalarField& baralpf_,
		const volVectorField& tildeUf_,
		const volVectorField& tildeUp_,		
		const volScalarField& barPg_,
		const volScalarField& p_,
		const volScalarField& rho_,
		const volVectorField& DragForce_,
		volVectorField& ResolvedDragForce_,
		vectorField& gii_
	);			
}	

#endif