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

\*---------------------------------------------------------------------------*/

#include "dsmcParcel.H"
#include "DsmcCloud.H"
#include "NoBinaryCollision.H"
#include "VariableHardSphere.H"
#include "LarsenBorgnakkeVariableHardSphere.H"
#include "ElasticChargedHardSphere.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef DsmcCloud<dsmcParcel> CloudType;

    makeBinaryCollisionModel(DsmcCloud<dsmcParcel>);

    // Add instances of collision model to the table
    makeBinaryCollisionModelType(NoBinaryCollision, CloudType);
    makeBinaryCollisionModelType(VariableHardSphere, CloudType);
    makeBinaryCollisionModelType(LarsenBorgnakkeVariableHardSphere, CloudType);
    makeBinaryCollisionModelType(ElasticChargedHardSphere, CloudType);

}


// ************************************************************************* //
