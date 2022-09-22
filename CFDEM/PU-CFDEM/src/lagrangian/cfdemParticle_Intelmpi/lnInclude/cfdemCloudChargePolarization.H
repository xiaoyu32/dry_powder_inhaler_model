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

#ifndef CfdemCloudChargePolarization_H
#define cfdemCloudChargePolarization_H

#include "cfdemCloud.H"
#include "wallDataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class chargeDensityModel;
class PolarizationModel;
class PermittivityModel;

/*---------------------------------------------------------------------------*\
                           Class cfdemCloudCharge Declaration
\*---------------------------------------------------------------------------*/

class CfdemCloudChargePolarization 
: 
public cfdemCloud
{
protected:
    
    autoPtr<chargeDensityModel> chargeDensityModel;
    autoPtr<PolarizationModel> polarizationModel;
    autoPtr<PermittivityModel> permittivityModel;
    
    // -- -- insulating walls -- --
    
    // -- cfd-dem wall data transfer model --
    autoPtr<wallDataExchangeModel> wallDataExchangeModel_;
    
    // -- dem wall surface charge density --
    mutable WallValueContainer wallSurfaceChargeDensity_;
    
    // -- cfd wall surface charge density --
    mutable volScalarField surfaceChargeDensity_;
    
    // -- normal electric field at surfaces (used for insulating wall charge transfer) --
    mutable WallValueContainer normalEfieldContainer_;
        
    // -- flag for insulating walls (track surface charge densities) --
    const bool insulatingWallsFlag_; 
     
    mutable double ** pcharge_;
    mutable double ** pef_;

    mutable double ** ef_grad;

public:

    // Constructors

    //- Construct from components
    CfdemCloudChargePolarization
    (
        const fvMesh& mesh
    );

    // Destructor

    ~CfdemCloudChargePolarization();

    
    class chargeDensityModel& chargeDensityM();
    class PolarizationModel& polarizationM();
    class PermittivityModel& permittivityM();
    
    // Member Functions
    void getDEMdata();    
    
    
    bool reAllocArrays() const;
	
    //bool evolve(volScalarField&);
    void evolveElectricField( volVectorField& Ef, volTensorField& gradEf );

    void evolve(
    		  volVectorField&,
                  volScalarField&,
                  volScalarField&,
                  volVectorField&,
                  volVectorField&,
		  volScalarField&
		);	

    inline wallDataExchangeModel& wallDataExchangeM()
    {
       return wallDataExchangeModel_();
    }
    
    inline bool insulatingWalls() const
    {
        return insulatingWallsFlag_;
    }
    	
    inline double ** charges() const
    {
        return pcharge_;
    };
	
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
