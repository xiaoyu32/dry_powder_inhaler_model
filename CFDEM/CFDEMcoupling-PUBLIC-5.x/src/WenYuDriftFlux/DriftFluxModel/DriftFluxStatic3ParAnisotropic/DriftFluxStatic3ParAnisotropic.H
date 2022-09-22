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
    Foam::dragModels::driftFluxModels::DriftFluxModel::DriftFluxStatic3ParAnisotropic

Description

SourceFiles
    DriftFluxStatic3ParAnisotropic.C

\*---------------------------------------------------------------------------*/

#ifndef DriftFluxStatic3ParAnisotropic_H
#define DriftFluxStatic3ParAnisotropic_H

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
                         Class DriftFluxStatic3ParAnisotropic Declaration
\*---------------------------------------------------------------------------*/
class DriftFluxStatic3ParAnisotropic
:
    public DriftFluxModel
{
  
    
  const dimensionedScalar filterSize_;
  const scalar lStar;
  const scalar uTerminal;
  const scalar alphaMax;
  const dimensionedVector g_;
  const scalar rhoParticle_;
  const scalar rhoFluid_;
  //mutable double DFFeatures;
  const int DFNFeatures;

  //- Reference to the mesh
  const fvMesh& mesh_;
  
  scalar uterminal_;
  
  //- Reference to carrier and solid phase properties
  mutable volScalarField solidVolumeFraction;
  mutable volVectorField uSlip;
  mutable volVectorField gradP;
  mutable volVectorField gradAlphaP;
  mutable volVectorField uSlipVdrift;
 
  //Neural network model
  mutable keras::KerasModel DFXnnModel_;
  mutable keras::KerasModel DFYnnModel_;
  mutable keras::KerasModel DFZnnModel_;

	
  mutable double* DFXmeans_;
  mutable double* DFXvars_;
  mutable double* DFYmeans_;
  mutable double* DFYvars_;
  mutable double* DFZmeans_;
  mutable double* DFZvars_;

   
public:

  //-Runtime type information
  TypeName("NeuralNet_Static3ParAnisotropic");
  
  //Constructor
  DriftFluxStatic3ParAnisotropic
  (
   const dictionary& dict,
   const phasePair& pair,
   const bool registerObject
   );

  // - Destructor
  virtual ~DriftFluxStatic3ParAnisotropic();
  
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
