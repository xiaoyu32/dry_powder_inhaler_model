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

Class
    forceModel

SourceFiles
    forceModel.C

\*---------------------------------------------------------------------------*/

#ifndef forceModel_H
#define forceModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

#include "probeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class forceModel Declaration
\*---------------------------------------------------------------------------*/

class forceModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        Switch treatExplicit_;                      // marker to treat force in implicit way (otherwise explicit)

        Switch treatDEM_;                           // marker to use the force only on DEM side

        Switch implDEM_;                           // marker to use the implicit force on DEM side

        mutable volVectorField impParticleForces_;  // sum of implicit particle forces [N]

        mutable volVectorField expParticleForces_;  // sum of explicit particle forces [N]

        bool coupleForce_;

        const word modelType_;

        mutable scalar cg_;
	
//_AO
        bool probeIt_;
//_AO
	
public:

    //- Runtime type information
    TypeName("forceModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            forceModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        forceModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~forceModel();


    // Selector

        static autoPtr<forceModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word forceType
        );


    // Member Functions
        virtual void setForce() const = 0;
	
	// -- communicate force model parameters to DEM --
	virtual void commValues() const
	{}
	
	// -- fetch data from DEM for force model --
	virtual void getDEMdata() const
	{}


	// -- allocate data for DEM exchange --
	virtual void reAllocArrays() const
	{}	
        //tmp<volScalarField> provideScalarField();

        virtual void manipulateScalarField(volScalarField&) const;

    // Access
        word modelType(){ return modelType_; };

        inline volVectorField& impParticleForces() const { return impParticleForces_;};

        inline volVectorField& expParticleForces() const { return expParticleForces_;};

        inline double ** impForces() const { return particleCloud_.impForces_;};

        inline double ** expForces() const { return particleCloud_.expForces_;};

        inline double ** DEMForces() const { return particleCloud_.DEMForces_;};

        inline double ** Cds() const { return particleCloud_.Cds_;};

        inline double ** fluidVel() const { return particleCloud_.fluidVel_;};

        inline const bool& coupleForce() const { return coupleForce_;};

//_AO 06/02/2014 
/*
        inline const scalar cg() const { return cg_; };

        void const setCG(double cg) const
        { 
            cg_ = cg;
            Info << "cg is set to: " << cg_ << endl;
        };
*/
//_AO
        void const setImpDEMdrag() const {particleCloud_.impDEMdrag_=true;};


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
