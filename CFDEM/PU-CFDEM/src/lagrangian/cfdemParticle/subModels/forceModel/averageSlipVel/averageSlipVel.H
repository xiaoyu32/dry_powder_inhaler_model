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

    calc vol average slip velocity and write to a field (postProc)

Class
    averageSlipVel

SourceFiles
    averageSlipVel.C

\*---------------------------------------------------------------------------*/

#ifndef averageSlipVel_H
#define averageSlipVel_H

#include "forceModel.H"
#include "momCoupleModel.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class averageSlipVel Declaration
\*---------------------------------------------------------------------------*/

class averageSlipVel
:
    public forceModel
{
private:
    dictionary propsDict_;

    word fluidVelFieldName_;        //fluid flow info
    const volVectorField& U_;

    word particleVelFieldName_;     //particle flow info
    const volVectorField& Us_;
    
    word  voidfractionFieldName_;   //void fraction field
    const volScalarField& voidfraction_;
    
    scalar rhoParticle_;
    word   rhoFluidName_;
    const  volScalarField& rhoFluid_;
    
    word   outputDirName_;
   
    OFstream* sPtrUSlip_;
    OFstream* sPtrIntMomentum_;
    OFstream* sPtrVelStats_;

public:

    //- Runtime type information
    TypeName("averageSlipVel");


    // Constructors

        //- Construct from components
        averageSlipVel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~averageSlipVel();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
