/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

Class
    Foam::dragModels::driftFluxModels::DriftFluxModel::DriftFluxNNetLocalNewX12Y6NoWCNoCap

Description

SourceFiles
    DriftFluxNNetLocalNewX12Y6NoWCNoCap.C

\*---------------------------------------------------------------------------*/

#ifndef DriftFluxNNetLocalNewX12Y6NoWCNoCap_H
#define DriftFluxNNetLocalNewX12Y6NoWCNoCap_H

#include "volFields.H"
#include "DriftFluxModel.H"
#include "dragModel.H"
#include "keras_model.H"
#include "extendedCentredCellToCellStencil.H"
#include "centredCPCCellToCellStencilObject.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
namespace driftFluxModels
{
  
  //#define XNFeatures 21
  //#define YNFeatures 12


/*---------------------------------------------------------------------------*\
                         Class DriftFluxNNetLocalNewX12Y6NoWCNoCap Declaration
\*---------------------------------------------------------------------------*/
class DriftFluxNNetLocalNewX12Y6NoWCNoCap
:
    public DriftFluxModel
{
  
  mutable volVectorField Ysgs;
  mutable volVectorField Xsgs;
  mutable volVectorField diffusion;
  mutable volVectorField sumXY;
  
  mutable volScalarField gradYsgs;  
  mutable volScalarField limiterVd;
    
  const Switch equilibriumSolutionSw_;  
    
    
  const dimensionedScalar driftFluxDiffusivity;
  const dimensionedScalar filterSize_;
  const scalar Xmin;
  const scalar Xmax;
  const scalar Ymin;
  const scalar Ymax;
  const scalar sumXYmax;
  const scalar sumXYmin;
  const scalar Vdmax;
  const scalar Vdmin;
  const scalar lStar;
  const scalar uTerminal;
  const scalar alphaMax;
  const dimensionedVector g_;
  const scalar rhoParticle_;
  const scalar rhoFluid_;
  //mutable double XNFeatures;
  //mutable double YNFeatures;
  const int XNFeatures;
  //const int XBNFeatures;
  const int YNFeatures;

  //- Reference to the mesh
  const fvMesh& mesh_;
  
  scalar uterminal_;
  
  //- Reference to carrier and solid phase properties
  mutable volScalarField solidVolumeFraction;
  mutable volVectorField uSlip;
  mutable volVectorField gradP;
  mutable volVectorField gradAlphaP;
  mutable volVectorField uSlipVdrift;
  mutable volVectorField driftFlux_min;

  //- Least-squares gradient vectors
  mutable List<List<scalar> > alphaList_;
  mutable List<List<vector> > uSlipList_;
  mutable List<List<vector> > gradPList_;
  mutable List<List<vector> > gradAlphaPList_;
  mutable List<List<vector> > driftFluxList_;
  //  mutable List<List<vector> > uSlipVdriftList_;
 

  //Neural network model
  mutable keras::KerasModel XnnModel_;
  mutable keras::KerasModel YnnModel_;
  //mutable keras::KerasModel XBnnModel_;
	
  mutable double* Xmeans_;
  mutable double* Xvars_;
  //mutable double* XBmeans_;
  //mutable double* XBvars_;
  mutable double* Ymeans_;
  mutable double* Yvars_;
   
public:

  //-Runtime type information
  TypeName("NeuralNet_LocalNew_X12Y6_NoWallCorrection_NoCap");
  
  //Constructor
  DriftFluxNNetLocalNewX12Y6NoWCNoCap
  (
   const dictionary& dict,
   const phasePair& pair,
   const bool registerObject
   );

  // - Destructor
  virtual ~DriftFluxNNetLocalNewX12Y6NoWCNoCap();
  
  virtual void evolve();

// Member Functions
	
      const fvMesh& mesh()
	{
	    return mesh_;
	}
	
	//- Return reference to the stencil
        const extendedCentredCellToCellStencil& stencil() const
        {
            return centredCPCCellToCellStencilObject::New(solidVolumeFraction.mesh());
	}
     
};

}
}
}
#endif
