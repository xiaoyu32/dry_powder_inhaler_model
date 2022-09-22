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

    -grad(p)*Vp * rho  ... if p is normalized with pressure
    -grad(p)*Vp        ... if p is real pressure
    pressure gradient force
    including interpolation of the velocity to the exact position

Class
    periodicPressure

SourceFiles
    periodicPressure.C

\*---------------------------------------------------------------------------*/

#ifndef periodicPressure_H
#define periodicPressure_H

#include "forceModel.H"
#include "momCoupleModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class periodicPressure Declaration
\*---------------------------------------------------------------------------*/

class periodicPressure
:
    public forceModel
{
private:
    dictionary propsDict_;

    bool verbose_;

    word velocityFieldName_;

    const volVectorField& U_;

    word densityFieldName_;

    const volScalarField& rho_;
       
    word voidfractionFieldName_;
        
    const volScalarField& voidfraction_;   

    word gravityFieldName_;  
          
    #if defined(version21) || defined(version16ext)
        const uniformDimensionedVectorField& g_;    // ref to gravity
    #elif defined(version15)
        const dimensionedVector& g_;    // ref to gravity
    #endif 
     	
    scalar rhoParticle_;
    
    scalar pressureMultiplier_;
    
public:

    //- Runtime type information
    TypeName("periodicPressure");


    // Constructors

        //- Construct from components
        periodicPressure
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~periodicPressure();


    // Member Functions
        void setForce() const;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
