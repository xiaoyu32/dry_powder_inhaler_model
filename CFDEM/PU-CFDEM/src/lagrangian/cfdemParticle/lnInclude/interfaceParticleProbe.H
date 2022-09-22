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

    calc interface velocity and write to a field (postProc)

Class
    interfaceParticleProbe

SourceFiles
    interfaceParticleProbe.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceParticleProbe_H
#define interfaceParticleProbe_H

#include "probeModel.H"
#include "fvCFD.H"
#include "polyMesh.H"
#include "cfdemCloud.H"
#include "OFstream.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceParticleProbel Declaration
\*---------------------------------------------------------------------------*/

class interfaceParticleProbe
:
    public probeModel
{
private:
    dictionary propsDict_;

    bool verbose_;

    bool twoDimensional_;

    bool depth_;

    word velFieldName_;

    const volVectorField& U_;

    word densityFieldName_;

    const volScalarField& rho_;

    word pressureFieldName_;

    const volScalarField& p_; 
   
    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    bool interpolation_; // use interpolated U field values
    
   // word   outputDirName_;
   
    //OFstream* sPtrInterfaceCheck_;
    OFstream* InterfaceCheckProps_;

public:

    //- Runtime type information
    TypeName("interfaceParticleProbe");


    // Constructors

        //- Construct from components
        interfaceParticleProbe
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word   typeName,
            char*  logFileName
        );

    // Destructor

        ~interfaceParticleProbe();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
