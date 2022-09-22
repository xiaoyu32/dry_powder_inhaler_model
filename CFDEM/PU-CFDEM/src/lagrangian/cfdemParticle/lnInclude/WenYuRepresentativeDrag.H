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

    Koch, Hill drag law
    based on Koch Hill 2001,"Inertial effects in suspensions and porous-media 
    flows", Annual Review of fluid mechanics.
    including interpolation of the velocity to the exact position
    including drag coefficient for implicit drag for DEM

Class
    WenYuDrag

SourceFiles
    WenYuDrag.C

\*---------------------------------------------------------------------------*/

#ifndef WenYuRepresentativeDrag_H
#define WenYuRepresentativeDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WenYuDrag Declaration
\*---------------------------------------------------------------------------*/

class WenYuRepresentativeDrag
:
    public forceModel
    
{
private:

    dictionary propsDict_;
    
    bool verbose_;
    
    word densityFieldName_; 
    
    word velFieldName_;

    word granVelFieldName_;

    word voidfractionFieldName_;
    
    word gravityFieldName_;
    
    const volScalarField& rho_;

    const volVectorField& U_;    // ref to fluid velocity

    const volVectorField& Us_;   // ref to granular velocity

    const volScalarField& voidfraction_;   // ref to voidfraction    

    bool interpolation_; // use interpolated field values
    
    bool interpolation_void_; // use interpolated voidage values only
    
    bool interpolation_vel_; // use interpolated velocity values only
       
// Parcel approach    
    bool parcelApproach_; 
    
    int npart_;
	
    scalar dpart_;
	
    scalar rhopart_;
    
    const int apiType_;
    
// Relative velocity at particle position: int<U-Us>@x_p)

    bool velslipdiff_;    

    bool weighting_; // use weighting for interpolation

    	
public:

    //- Runtime type information
    TypeName("WenYuRepresentativeDrag");


    // Constructors

        //- Construct from components
        WenYuRepresentativeDrag
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~WenYuRepresentativeDrag();


    // Member Functions
        void setForce() const;
	
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
