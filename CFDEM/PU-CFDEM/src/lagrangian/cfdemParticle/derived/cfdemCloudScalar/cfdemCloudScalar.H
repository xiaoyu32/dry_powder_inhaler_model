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

#ifndef CfdemCloudScalar_H
#define CfdemCloudScalar_H

#include "cfdemCloud.H"
#include "passiveScalarDensityModel.H"
#include "List.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class PassiveScalarDensityModel;

/*---------------------------------------------------------------------------*\
                           Class CfdemCloudScalar Declaration
\*---------------------------------------------------------------------------*/

class CfdemCloudScalar 
: 
public cfdemCloud
{
protected:
    
    autoPtr<passiveScalarDensityModel> passiveScalarDensityModel_;
    
    
    //change of passive scalar on particle
    mutable List<double **> dpscalar_;
    
    //passive scalar on particles
    mutable List<double **> pscalar_;
      
    // Velocity correction for periodic domain
    mutable double ** velCorrection_;
    
    // Collisional stresses
    mutable double ** collStress_;
    volSymmTensorField sigmaCollCell_;  
    
    mutable int numberOfPassiveScalars_;
    
public:

    // Constructors

    //- Construct from components
    CfdemCloudScalar
    (
        const fvMesh& mesh
    );

    // Destructor

    ~CfdemCloudScalar();

    
    passiveScalarDensityModel& passiveScalarDensityM();

    // Member Functions
    void getDEMdata();

    bool reAllocArrays() const;
	

    void evolvePassiveScalar(	
    				List< autoPtr<volScalarField> >&,
                      		List< autoPtr<volScalarField> >&,
                      		List< autoPtr<volScalarField> >&,
				volScalarField&,
                      		volVectorField&,
                      		volVectorField&
			    );	
    
    void updatePassiveScalarSource( List< autoPtr<volScalarField> >& );
    
	
    void calcSendVelocityCorrection
    (
	scalar&,
	volScalarField&,
	volScalarField&,
	volVectorField&,
	volVectorField&
    );	

    void getDEMStressData();

    symmTensor collStress(int);

    void writeCollStress
    (
	fvMesh&,
	cfdemCloud&	
    );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
