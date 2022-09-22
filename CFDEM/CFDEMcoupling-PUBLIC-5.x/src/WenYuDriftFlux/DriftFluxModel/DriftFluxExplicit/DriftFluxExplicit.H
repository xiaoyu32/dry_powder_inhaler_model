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
    Foam::dragModels::driftFluxModels::DriftFluxModel::DriftFluxExplicit

Description

SourceFiles
    DriftFluxExplicit.C

\*---------------------------------------------------------------------------*/

#ifndef DriftFluxExplicit_H
#define DriftFluxExplicit_H

#include "volFields.H"
#include "DriftFluxModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
namespace driftFluxModels
{

/*---------------------------------------------------------------------------*\
                         Class DriftFluxExplicit Declaration
\*---------------------------------------------------------------------------*/
class DriftFluxExplicit
:
    public DriftFluxModel
{
    
  const dimensionedScalar driftFluxDiffusivity;
  const dimensionedScalar filterSize_;
  const dimensionedVector g_;
  
  const scalar XsgsCoeff; 
  const scalar YsgsCoeff;
    
public:

  //-Runtime type information
  TypeName("Explicit");
  
  //Constructor
  DriftFluxExplicit
  (
   const dictionary& dict,
   const phasePair& pair,
   const bool registerObject
   );

  // - Destructor
  virtual ~DriftFluxExplicit();
  
  virtual void evolve();
  
};

}
}
}
#endif
