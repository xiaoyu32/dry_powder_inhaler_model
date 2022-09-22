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
    Foam::dragModels::driftFluxModels::DriftFluxModel::OzelIJMF

Description

SourceFiles
    OzelIJMF.C

\*---------------------------------------------------------------------------*/

#ifndef OzelIJMF_H
#define OzelIJMF_H

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
                         Class OzelIJMF Declaration
\*---------------------------------------------------------------------------*/
class OzelIJMF
:
    public DriftFluxModel
{
    
public:

  //-Runtime type information
  TypeName("OzelIJMF");
  
  //Constructor
  OzelIJMF
  (
   const dictionary& dict,
   const phasePair& pair,
   const bool registerObject
   );

  // - Destructor
  virtual ~OzelIJMF();
  
  virtual void evolve();
  
};

}
}
}
#endif
