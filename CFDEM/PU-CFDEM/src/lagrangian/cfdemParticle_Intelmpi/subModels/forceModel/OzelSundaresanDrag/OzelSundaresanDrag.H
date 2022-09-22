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

    Ozel / Sundaresan drag law
    - based on the Wen-Yu drag model (with the option to switch to Ergun
      for dense situations)
    - allows the use of average particle velocity to be used 
      in the calculation of the slip velocity between particles and fluid
    - includes fluid-size drag correction
    - allows the specification of an-isotropic drag coefficient
    - including interpolation of the velocity to the exact position

Copyright of this Contribution
    Copyright 2014-     Ali Ozel & Sankaran Sundaresan, Princeton University

Class
    OzelSundaresanDrag

SourceFiles
    OzelSundaresanDrag.C

\*---------------------------------------------------------------------------*/

#ifndef OzelSundaresanDrag_H
#define OzelSundaresanDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class OzelSundaresanDrag Declaration
\*---------------------------------------------------------------------------*/

class OzelSundaresanDrag
:
    public forceModel
{
private:
    dictionary propsDict_;

    bool verbose_;

    word velFieldName_;

    const volVectorField& U_;

    word granVelFieldName_;

    const volVectorField& Us_;   // ref to granular velocity

    word densityFieldName_;

    const volScalarField& rho_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

//_AO 06/02/2014 
//    const scalar phi_;
//_AO 06/02/2014 
  
    bool interpolation_; // use interpolated field values

    mutable scalar scaleDia_;

    mutable scalar scaleDrag_;

    mutable scalar switchingVoidfraction_; //voidfraction above which dilute formulation will be used

// Use Relative velocity at particle position: int<U-Us>@x_p)

    bool velslipdiff_;  
    
// Wen & Yu with Ergun's correlation

    bool ErgunMix_;    
    
// Model constants
    mutable  scalar Kxx_;
    mutable  scalar Kyy_;
    mutable  scalar Kzz_;    
    mutable  scalar const_h1;  
    mutable  scalar const_h2; 
    mutable  scalar const_h3;
    mutable  scalar alpsmax;  

    word gravityFieldName_;

    #ifdef version21
        const uniformDimensionedVectorField& g_;    // ref to gravity
    #elif defined(version16ext) || defined(version15)
        const dimensionedVector& g_;    // ref to gravity
    #endif

public:

    //- Runtime type information
    TypeName("OzelSundaresanDrag");


    // Constructors

        //- Construct from components
        OzelSundaresanDrag
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~OzelSundaresanDrag();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
