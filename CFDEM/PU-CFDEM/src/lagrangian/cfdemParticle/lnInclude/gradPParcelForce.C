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
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "gradPParcelForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradPParcelForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    gradPParcelForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
gradPParcelForce::gradPParcelForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    pFieldName_(propsDict_.lookup("pFieldName")),
    p_(sm.mesh().lookupObject<volScalarField> (pFieldName_)),
    velocityFieldName_(propsDict_.lookup("velocityFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velocityFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    useRho_(false),
    useU_(false),
    interpolation_(false),
    weighting_(false),
    parcelApproach_(false),
    npart_(1),
    dpart_(1.0),
    rhopart_(1.0),
    apiType_
    (
        propsDict_.lookupOrDefault("apiType", 1)
    )   
{
    if (modelType_ == "B")
    {
        FatalError <<"using  model gradPForce with model type B is not valid\n" << abort(FatalError);
    }else
    {
        treatDEM_=true;
        Info << "gradPForce is applied only to DEM side" << endl;
    }
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    if (propsDict_.found("useU")) useU_=true;
    if (propsDict_.found("interpolation"))
    {
        Info << "using interpolated value of pressure gradient." << endl;
        interpolation_=true;
    }

    if(p_.dimensions()==dimensionSet(0,2,-2,0,0))
        useRho_ = true;

    particleCloud_.checkCG(true);
    
    if (propsDict_.found("weighting")) 
    {
    	weighting_= true;    
	interpolation_= false;
        Info << " Weigthing is using for interpolation " << endl;
    }
    
    // Parcel Approach
    parcelApproach_=true;
    npart_  = readScalar(propsDict_.lookup("npart"));
    dpart_  = readScalar(propsDict_.lookup("dpart"));  
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gradPParcelForce::~gradPParcelForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gradPParcelForce::setForce() const
{
    volVectorField gradPField = fvc::grad(p_);
    /*if (useU_)
    {
        // const volScalarField& voidfraction_ = particleCloud_.mesh().lookupObject<volScalarField> ("voidfraction");
        volScalarField U2 = U_&U_;// *voidfraction_*voidfraction_;
        if (useRho_)
            gradPField = fvc::grad(0.5*U2);
        else
            gradPField = fvc::grad(0.5*rho_*U2);
    }*/
//
//
    vector gradP;
    scalar ds;
    scalar Vs;
    scalar rho;
    vector position;
    vector force;
    label cellI;

    interpolationCellPoint<vector> gradPInterpolator_(gradPField);


    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            force=vector(0,0,0);
            cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                position = particleCloud_.position(index);

                if(interpolation_) // use intepolated values for alpha (normally off!!!)
                {
                    gradP = gradPInterpolator_.interpolate(position,cellI);
                }
		else if(weighting_)
		{
		    for (int subCell = 0; subCell < particleCloud_.voidFractionM().cellsPerParticle()[index][0]; subCell++)	
		    {
		    	gradP += particleCloud_.weighting(index,subCell)*gradPField[cellI];
		    } 	
		}else
                {
                    gradP = gradPField[cellI];
                }

                ds = 2*particleCloud_.radius(index);
		
		// Parcel approach: based on particle diameter
		int ptype = particleCloud_.types()[index][0];
		
		if ( (parcelApproach_) && (ptype == apiType_) )
		{
	 	    ds = dpart_;
		}
		
                Vs = ds*ds*ds*M_PI/6;
                rho = rho_[cellI];

                // calc particle's pressure gradient force
                if (useRho_)
                    force = -Vs*gradP*rho;
                else
                    force = -Vs*gradP;
		    
		// Parcel approach: drag on one particle multiply by number of particles in a parcel
		if ( (parcelApproach_) && (ptype == apiType_) )
		{
		    force *= npart_;
		}

                if(verbose_ && index >=0 && index <1)
                {
                    //Info << "periodicPressure " << endl; 
		    //Info << "index = " << index << endl;
                    Info << "gradP = " << gradP << endl;
                    Info << "force = " << force << endl;
                }
            }

            // set force on particle
            if(!treatDEM_){
                if(!treatExplicit_) for(int j=0;j<3;j++) impForces()[index][j] += force[j];
                else  for(int j=0;j<3;j++) expForces()[index][j] += force[j];
            }
            for(int j=0;j<3;j++) DEMForces()[index][j] += force[j];

        //}
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
