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

#include "GaussVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(GaussVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    GaussVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
GaussVoidFraction::GaussVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    alphaLimited_(0),
    scaleUpVol_(readScalar(propsDict_.lookup("scaleUpVol"))),
    verbose_(false)
{
    //Info << "\n\n W A R N I N G - do not use in combination with differentialRegion model! \n\n" << endl;
    //Info << "\n\n W A R N I N G - this model does not yet work properly! \n\n" << endl;
    //reading maxCellsPerParticle from dictionary
    maxCellsPerParticle_= 27; //readLabel(propsDict_.lookup("maxCellsPerParticle"));

    if(scaleUpVol_ < 1){ FatalError<< "scaleUpVol shloud be > 1."<< abort(FatalError); }
    if(alphaMin_ > 1 || alphaMin_ < 0.01){ FatalError<< "alphaMin shloud be > 1 and < 0.01." << abort(FatalError); }
    if (propsDict_.found("weight"))
        setWeight(readScalar(propsDict_.lookup("weight")));
    if (propsDict_.found("verbose")) verbose_ = true;
        
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

GaussVoidFraction::~GaussVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void GaussVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes) const
{
    reAllocArrays();

    voidfractionNext_.internalField() = 1;
    
    // Neighbouring cells
    CPCCellToCellStencil neighbourCells(particleCloud_.mesh());    

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            //reset
            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell]=0;
                particleVolumes[index][subcell]=0;
            }
            cellsPerParticle_[index][0]=1.0;

            //collecting data
            scalar pi = M_PI;
            label particleCenterCellID=particleCloud_.cellIDs()[index][0];
            scalar radius =  particleCloud_.radii()[index][0];
              //Info << "physical radius" << radius  << endl;

            scalar volume = 4./3.*radius*radius*radius*pi*weight();
            radius = radius*pow(scaleUpVol_,0.3333);
              //Info << "fictuous radius" << radius  << endl;
            vector positionCenter=particleCloud_.position(index);
	    scalar core;
	    scalar dist;

	    scalarField particleWeightsGauss(maxCellsPerParticle_,scalar(0.0));
	    scalar dist_s(0);
	    scalar sumWeights(0);
	    
	    // debug
            if(verbose_ && index==0) Info << "Particle 0 cell ID  = " << particleCenterCellID << endl; 

            if (particleCenterCellID >= 0)
            {
                //labelHashSet hashSett;

                //determining label and degree of coveredness of cells covered by the particle
                //buildLabelHashSet(radius*3.0, positionCenter, particleCenterCellID, hashSett);
                  //Info << "completeSize=" << hashSett.size() << ", completeList =\n" << endl;
                  //for(label i=0;i<hashSett.size();i++) Info << " ," << hashSett.toc()[i] << endl;	
                 
		/* 
                //generating list with cell and subcells
                scalar hashSetLength = hashSett.size();
                if (hashSetLength > maxCellsPerParticle_)
                {
                    FatalError<< "big particle algo found more cells ("<< hashSetLength
                              <<") than storage is prepared ("<<maxCellsPerParticle_<<")" << abort(FatalError);
                }
                else if (hashSetLength > 0)
                {
                    cellsPerParticle_[index][0]=hashSetLength;

                    //making sure that the cell containing the center is the first subcell
                    particleCloud_.cellIDs()[index][0]=particleCenterCellID;
                    //deleting the cell containing the center of the particle
                    hashSett.erase(particleCenterCellID);

                    //==========================//
                    //setting the voidfractions

		    dist = mag(particleCloud_.mesh().C()[particleCenterCellID]-particleCloud_.position(index));
		    core = pow(2.0/radius/radius/M_PI,1.5)*exp(-dist*dist/2.0/radius/radius)*particleCloud_.mesh().V()[particleCenterCellID];

                    // volume occupied in every covered cell
                    scalar occupiedVolume = volume*core;

                    // correct volumefraction of centre
                    voidfractionNext_[particleCenterCellID] -=occupiedVolume/particleCloud_.mesh().V()[particleCenterCellID];

                    particleWeights[index][0] += core;
                    particleVolumes[index][0] += occupiedVolume;

                      //Info << "Centre:set voidfraction in cellI=" << particleCenterCellID
                      //     << ", voidfraction =" << voidfractionNext_[particleCenterCellID] << endl;

                    // correct volumefraction of sub-cells
                    for(label i=0;i<hashSetLength-1;i++)
                    {
                        label cellI=hashSett.toc()[i];
                        particleCloud_.cellIDs()[index][i+1]=cellI; //adding subcell represenation

		        dist = mag(particleCloud_.mesh().C()[cellI]-particleCloud_.position(index));
		        core = pow(2.0/radius/radius/M_PI,1.5)*exp(-dist*dist/2.0/radius/radius)*particleCloud_.mesh().V()[cellI];
                        scalar occupiedVolume = volume*core;
                        voidfractionNext_[cellI] -=occupiedVolume/particleCloud_.mesh().V()[cellI];
                        particleWeights[index][i+1] += core;
                        particleVolumes[index][i+1] += occupiedVolume;

                          //Info << "AFTER:set voidfraction in cellI=" << cellI
                          //     << ", voidfraction =" << voidfractionNext_[cellI] << endl;

                    }

                    // debug
                    if(index==0)
                    {
                        Info << "particle 0 is represented by " << hashSetLength << "cells" << endl;
                    }
		    */
		    
		    const labelList& cellsNeigh = particleCloud_.mesh().cellCells()[particleCenterCellID]; //neighbourCells[particleCenterCellID];
		    sumWeights = 0;
		    dist_s = 0;
		    
		    cellsPerParticle_[index][0] = maxCellsPerParticle_;
		    
		    //making sure that the cell containing the center is the first subcell
                    particleCloud_.cellIDs()[index][0] = particleCenterCellID;

		    //Info << " index = " << index << " ii = " << ii << " cellID = " << cellID << endl;

		    forAll(cellsNeigh,jj)
		    {
			    if(jj!=0) particleCloud_.cellIDs()[index][jj] = cellsNeigh[jj]; //adding subcell represenation
						    
			    // Find distances between particle and neighbouring cells					
			    dist_s = mag(particleCloud_.mesh().C()[cellsNeigh[jj]]-positionCenter)/pow(particleCloud_.mesh().V()[cellsNeigh[jj]],1./3.);

			    if(dist_s <= 0.5)
			    {		
				    particleWeightsGauss[jj] =  1./4.*pow(dist_s,4)-5./8.*pow(dist_s,2)+115./192.;
			    }
			    else if (dist_s > 0.5 && dist_s <= 1.5)
			    {		
				    particleWeightsGauss[jj] = -1./6.*pow(dist_s,4)+5./6.*pow(dist_s,3)-5./4.*pow(dist_s,2)+5./24.*dist_s+55./96.;
			    }
			    else if (dist_s > 1.5 && dist_s <= 2.5)
			    {		
				    particleWeightsGauss[jj] =  pow(2.5-dist_s,4)/24.;
			    }
			    else
			    {		
				    particleWeightsGauss[jj] = 0;
			    }

			    sumWeights += particleWeightsGauss[jj];
			    
			    // debug
                    	    if(verbose_ && index==0) Info << "Particle 0 gaussian weights = " << particleWeightsGauss[jj] << endl;

		    }
		    // debug
                    if(verbose_ && index==0) Info << "Particle 0 sum of gaussian weights = " << sumWeights << endl;
    
		    forAll(cellsNeigh,jj)
		    {	
			    if ( sumWeights != 0 )
			    {
                		    // volume occupied in every covered cell
                		    scalar occupiedVolume = volume*particleWeightsGauss[jj]/sumWeights;

                		    // correct volumefraction of centre
                		    voidfractionNext_[cellsNeigh[jj]] -= occupiedVolume/particleCloud_.mesh().V()[cellsNeigh[jj]];
				    
				    if(verbose_ && index==0)
                    		    {   
					Info << "Particle 0 void fraction contribution into cell " << cellsNeigh[jj] 
					     << " " << voidfractionNext_[cellsNeigh[jj]] << endl;
                    		    }	

                		    particleWeights[index][jj] = particleWeightsGauss[jj]/sumWeights;
                		    particleVolumes[index][jj] = occupiedVolume;
			    }
			    else
			    {
			    	    particleWeights[index][jj] = 0.;
                		    particleVolumes[index][jj] = 0.;	
			    }
		    }
							    
                    //==========================//
                //}//end cells found on this proc
            }// end found cells
        //}// end if masked
    }// end loop all particles

    //bringing eulerian field to particle array
    for(label index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        for(label subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
        {
            label cellID = particleCloud_.cellIDs()[index][subcell];

            if(cellID >= 0)
            {
                 // limiting voidfraction
                 if (voidfractionNext_[cellID] < alphaMin_) voidfractionNext_[cellID]=alphaMin_;

                 // set particle based voidfraction
                 voidfractions[index][subcell] = voidfractionNext_[cellID];
                 //Info<<"setting the voidfraction, index = "<<index<<endl;
            }
            else
            {
                voidfractions[index][subcell] = -1.;
            }
        }
    }
}

void GaussVoidFraction::buildLabelHashSet
(
    const scalar radius,
    const vector position,
    const label cellID,
    labelHashSet& hashSett
)const
{
    hashSett.insert(cellID);
    //Info<<"cell inserted"<<cellID<<endl;
    const labelList& nc = particleCloud_.mesh().cellCells()[cellID];
    forAll(nc,i){
        label neighbor=nc[i];
        if(!hashSett.found(neighbor) && mag(position-particleCloud_.mesh().C()[neighbor])<radius){
            buildLabelHashSet(radius,position,neighbor,hashSett);
        }
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
