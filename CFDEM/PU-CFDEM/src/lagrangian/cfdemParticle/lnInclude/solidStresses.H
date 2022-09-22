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

    solids pressure for MP-PIC, following Benyaha and Sundresan (2012),PowTech

Class
    solidStresses

SourceFiles
    solidStresses.C

\*---------------------------------------------------------------------------*/

#ifndef solidStresses_H
#define solidStresses_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidStresses Declaration
\*---------------------------------------------------------------------------*/

class solidStresses
:
    public forceModel
{
private:
    dictionary propsDict_;
    
    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    word solidVelocityFieldName_;
    
    const volVectorField& Us_;
    
    const scalar phic_;    

    bool interpolation_; 

    mutable volScalarField Pp_;

    mutable volVectorField gradPp_;

    bool onlyPp_;  

    bool blendedModel_;

    mutable volScalarField mup_;

    mutable volVectorField gradTaup_;
    
    bool verbose_;
    
    scalar rhoParticle_;
    
    scalar dp_;
    
    scalar phid_;

public:

    //- Runtime type information
    TypeName("solidStresses");


    // Constructors

        //- Construct from components
        solidStresses
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~solidStresses();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
