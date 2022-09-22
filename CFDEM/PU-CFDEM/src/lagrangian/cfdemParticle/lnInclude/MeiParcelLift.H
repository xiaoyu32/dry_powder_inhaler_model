/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

    Mei Lift force following Loth and Dorgan 2009
    contribution from RQ

Class
    MeiLift

SourceFiles
    MeiLift.C

\*---------------------------------------------------------------------------*/

#ifndef MeiParcelLift_H
#define MeiParcelLift_H

#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MeiLift Declaration
\*---------------------------------------------------------------------------*/

class MeiParcelLift
:
    public forceModel
{
private:
    dictionary propsDict_;

    bool verbose_;

    word velFieldName_;

    const volVectorField& U_;

    word densityFieldName_;

    const volScalarField& rho_;

    /*word vorticityFieldName_;

    volVectorField& vorticity_;*/
    
    // Parcel approach
    bool parcelApproach_;
    
    int npart_;
    
    scalar dpart_;
    
    scalar rhopart_;

    const int apiType_;

public:

    //- Runtime type information
    TypeName("MeiParcelLift");


    // Constructors

        //- Construct from components
        MeiParcelLift
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~MeiParcelLift();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
