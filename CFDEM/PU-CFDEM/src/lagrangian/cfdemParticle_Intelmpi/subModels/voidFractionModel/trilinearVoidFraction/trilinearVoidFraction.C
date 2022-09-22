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

#include "trilinearVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"

//#include "mpi.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(trilinearVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    trilinearVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
trilinearVoidFraction::trilinearVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    alphaLimited_(0),
    tooMuch_(0.0),
    scaleUpVol_(readScalar(propsDict_.lookup("scaleUpVol"))),
    interpolation_(false)
{
    maxCellsPerParticle_ = 8;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

trilinearVoidFraction::~trilinearVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void trilinearVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes) const
{
    reAllocArrays();

    scalar pi = M_PI;
    vector position(0,0,0);
    label  cellID=-1;
    scalar radius(-1);
    scalar Volp(0);
    scalar partCellVol(0);
 
    /*
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            // reset
            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell]=0;
                particleVolumes[index][subcell]=0;
            }

	    //cellsPerParticle()[index][0] = 1;
            cellsPerParticle()[index][0] = 0;
	    position = particleCloud_.position(index);
            cellID = particleCloud_.cellIDs()[index][0];
            radius = particleCloud_.radii()[index][0];
            Volp =  4./3.*radius*radius*radius*pi*weight();
            radius = radius*pow(scaleUpVol_,1/3);
            partCellVol = 0;
	    scalar newAlpha = 0;	

            if (cellID > -1)  // particel centre is in domain
            {
				
		label partCellId = cellID;
		
		// Cell to vertices
		labelList& CelltoPoints = particleCloud_.mesh().cellPoints();//[cellID];
		
		// Find the closest vertice
		scalar dist_max(1.e08);
		scalar dist(0);
		label pImin;
		forAll(CelltoPoints,pI)
		{
			label point = CelltoPoints[pI];
			dist = mag(position-particleCloud_.mesh().points()[point]);
			if(dist<dist_max) 
			{
				dist_max = dist;
				pImin = pI;
			}
		}				
		
		Info << " Closest vertice = " << pImin << endl;
		
		// From closest vertice to obtain supporting cells
		
		// int=0 the cell which particle belongs to ...
		for(int i=0; i < maxCellsPerParticle_ ; i++)
        	{	
			particleCloud_.cellIDs()[index][i] = partCellId;
			
			particleWeights[index][i] = 1/static_cast<scalar>(maxCellsPerParticle_) ; //particleCloud_.mesh().V()[cellI]
			particleVolumes[index][i] = Volp*particleWeights[index][i];
		
			partCellVol = particleCloud_.mesh().V()[partCellId];
        		newAlpha += particleVolumes[index][i] / partCellVol;	
			
			cellsPerParticle()[index][0]++;		
		}
		
        	if( (voidfractionNext_[cellID] - newAlpha) > alphaMin_) 
		{
			voidfractionNext_[cellID] -= newAlpha;
        	}
		else
        	{
        	    	voidfractionNext_[cellID] = alphaMin_;
        	    	tooMuch_ += (alphaMin_-newAlpha) * particleCloud_.mesh().V()[cellID];
        	}		
     */
     /*
                //OUTPUT
                if (index==0)
                {					
                    for(int i=0; i < cellsPerParticle()[index][0] ; i++)
                    {
			Info << index << " subCell = " << i << " CellI = " << particleCloud_.cellIDs()[index][i] 
				      << " Particle weighting = " << particleWeights[index][i]  << endl;
		    }
                }
		//
						
            }
	    
        if(index == particleCloud_.numberOfParticles()-1) Info << "Total particle volume neglected: " << tooMuch_<< endl;
    }// end loop all particles


    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
	for(int subcell=0 ; subcell < cellsPerParticle()[index][0] ; subcell++)
	{
		label cellID = particleCloud_.cellIDs()[index][subcell];

                if(cellID >= 0)
                {
                    voidfractions[index][subcell] = voidfractionNext_[cellID];
                } 
                else
                {
                    voidfractions[index][subcell] = -1.;
		}
	}
    }
    */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
