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

    Gaussian distributed particle voidfraction model
    contribution from RQ

Class
    GaussVoidFraction

SourceFiles
    GaussVoidFraction.C

\*---------------------------------------------------------------------------*/

#ifndef GaussVoidFraction_H
#define GaussVoidFraction_H

#include "voidFractionModel.H"
#include "CPCCellToCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class GaussVoidFraction
:
    public voidFractionModel
{

private:
        dictionary propsDict_;

        const scalar alphaMin_;          //NP min value of voidFraction

	mutable bool alphaLimited_;

        const scalar scaleUpVol_;        //NP scaling radius, keeping volume of particle
	
	mutable bool verbose_;

public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct from components
        GaussVoidFraction
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~GaussVoidFraction();


    // Member Functions
        void setvoidFraction(double** const& ,double**&, double**&, double**&) const;

        void buildLabelHashSet
        (
            const scalar radius,
            const vector position,
            const label cellID,
            labelHashSet& hashSett
            /*labelHashSet& hashSetComplete,
            labelHashSet& hashSetPartial,
            labelHashSet& hashSetRest*/
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
