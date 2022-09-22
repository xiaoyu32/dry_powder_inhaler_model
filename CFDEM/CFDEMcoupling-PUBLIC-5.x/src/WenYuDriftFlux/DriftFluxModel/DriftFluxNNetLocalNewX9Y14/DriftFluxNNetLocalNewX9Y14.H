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
    Foam::dragModels::driftFluxModels::DriftFluxModel::DriftFluxNNetLocalNewX9Y14

Description

SourceFiles
    DriftFluxNNetLocalNewX9Y14.C

\*---------------------------------------------------------------------------*/

#ifndef DriftFluxNNetLocalNewX9Y14_H
#define DriftFluxNNetLocalNewX9Y14_H

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
                         Class DriftFluxNNetLocalNewX9Y14 Declaration
\*---------------------------------------------------------------------------*/
class DriftFluxNNetLocalNewX9Y14
:
    public DriftFluxModel
{
  
  mutable volVectorField Ysgs;
  mutable volVectorField Xsgs;
  mutable volVectorField sumXY;
    
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
  const int YNFeatures;

  //- Reference to the mesh
  const fvMesh& mesh_;
  
  scalar uterminal_;
  
  //- Reference to carrier and solid phase properties
  /*mutable volScalarField solidVolumeFraction;
  mutable volVectorField uSlip;
  mutable volVectorField gradP;
  mutable volVectorField gradGradPz;
  mutable volVectorField gradUslipZ;
  mutable volVectorField gradAlphaP;
  mutable volVectorField uSlipVdrift;
  mutable volVectorField driftFlux_min;*/

  //- Least-squares gradient vectors
  /*mutable List<List<scalar> > alphaList_;
  mutable List<List<vector> > uSlipList_;
  mutable List<List<vector> > gradPList_;
  mutable List<List<vector> > gradAlphaPList_;
  mutable List<List<vector> > driftFluxList_;*/
  //  mutable List<List<vector> > uSlipVdriftList_;
 

  //Neural network model
  mutable keras::KerasModel XnnModel_;
  mutable keras::KerasModel YnnModel_;
	
  mutable double* Xmeans_;
  mutable double* Xvars_;
  mutable double* Ymeans_;
  mutable double* Yvars_;
   
public:

  //-Runtime type information
  TypeName("NeuralNet_LocalNew_X9Y14");
  
  //Constructor
  DriftFluxNNetLocalNewX9Y14
  (
   const dictionary& dict,
   const phasePair& pair,
   const bool registerObject
   );

  // - Destructor
  virtual ~DriftFluxNNetLocalNewX9Y14();
  
  virtual void evolve();

// Member Functions
	
      const fvMesh& mesh()
	{
	    return mesh_;
	}
     
};

}
}
}
#endif
