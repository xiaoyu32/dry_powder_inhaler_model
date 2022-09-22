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

    cloud class managing DEM data for CFD-DEM coupling and IB representation

Class
    Foam::cfdemCloudIB derived from cfdemCloud

SourceFiles
    cfdemCloudIB.C

\*---------------------------------------------------------------------------*/

#ifndef cfdemCoarseDPICloud_H
#define cfdemCoarseDPICloud_H

#include "cfdemCloud.H"
#include "wallDataExchangeModel.H"
#include "apiTransferModel.H"
#include "apiDensityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//class apiTransferModel;
//class apiDensityModel;
//class wallDataExchangeModel;
//class WallValueContainer;

/*---------------------------------------------------------------------------*\
                           Class cfdemCloudCharge Declaration
\*---------------------------------------------------------------------------*/

class cfdemCoarseDPICloud 
: 
public cfdemCloud
{
protected:
    
    // -- carrier particle API particle concentration --
    mutable double ** carrierConcentration_;
    mutable double ** carrierFlux_;
    
    // -- transfer model for computing API particle interphase transfer (no collisions) --
    autoPtr<apiTransferModel> apiTransferModel_;
    
    // -- api property interpolation model --
    autoPtr<apiDensityModel> apiDensityModel_;
    
    // -- wall data exchange model --
    autoPtr<wallDataExchangeModel> wallDataExchangeModel_;
    
    // -- references to continuous fields -- 
    const volScalarField& alpha_;
    const volVectorField& Us_;
    const volVectorField& U_;
    
    // -- dem wall concentration --
    mutable WallValueContainer demWallConcentration_;
    
    // -- continuous API particle fields --
    mutable volScalarField apiFluidDensity_;
    mutable volScalarField apiCarrierDensity_;
    
    // -- coefficient of fluid reatachment and detachement fluxes --
    mutable volScalarField apiImplSource_;
    
    // -- api particle source/sink due to carrier particle and wall collisions -- 
    mutable volScalarField apiExplSource_;
    
    
    // -- deposited api particle surface density at wall --
    mutable surfaceScalarField apiWallDensity_;
    
    // -- api particle dispersion coefficient --
    mutable volScalarField apiDiffusivity_;
    
    // -- evolve api density field --
    virtual void evolveAPIfields();
    
    virtual void updateExplSources();
    virtual void updateImplSources();
    
public:

    // Constructors

    //- Construct from components
    cfdemCoarseDPICloud
    (
        const fvMesh& mesh,
	const volScalarField& alpha_,
	const volVectorField& Us_,
	const volVectorField& U_
    );

    // Destructor

    ~cfdemCoarseDPICloud();
    
    inline surfaceScalarField& apiWallDensity()
    {
        return apiWallDensity_;
    }

    inline volScalarField& apiCarrierDensity()
    {
        return apiCarrierDensity_;
    }    
    
    inline apiDensityModel& apiDensityM()
    {
       return apiDensityModel_();
    }
    
    inline wallDataExchangeModel& wallDataExchangeM()
    {
       return wallDataExchangeModel_();
    }
 
    inline apiTransferModel& apiTransferM()
    {
       return apiTransferModel_();
    } 
    
    
    virtual void getDEMdata();
    virtual bool reAllocArrays() const;
	
    //bool evolve(volScalarField&);
    virtual void commValues();

    virtual bool evolve(
                 volScalarField&,
                 volVectorField&,
                 volVectorField&
	       );	

    inline double ** carrierFlux() const
    {
        return carrierFlux_;
    };
	
    inline double ** carrierConcentration() const
    {
        return carrierConcentration_;
    };
	
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
