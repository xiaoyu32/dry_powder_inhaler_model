{
    scalar pi = M_PI;
    vector position(0,0,0);
    label cellID(-1);
    scalar radius(-1.);
    scalar cellVol(0);
    scalar tooMuch(0.);
    
    // Lists for mapping function
    cellsPerParticle.setSize(nP[me]);
    neighboringCellIDs.setSize(nP[me],subCellIDsList);
    
    particleWeights.setSize(nP[me],subCellList);
    particleVolumes.setSize(nP[me],subCellList);
    
    if(!notWeighthing)
    {
	    for(label index = 0; index < nP[me]; index++)
	    {

        	    cellsPerParticle[index] = 1;
        	    position = positions[me][index];
        	    cellID = cellIDs[index];
        	    radius = radii[me][index];
        	    scalar volume =  4./3.*pi*radius*radius*radius;
        	    cellVol = 0.;

        	    //--variables for sub-search
        	    int nPoints = 29;
        	    int nNotFound=0,nUnEqual=0,nTotal=0;
        	    vector offset(0,0,0);
        	    int cellsSet = 0;

        	    if (cellID >-1)  // particel centre is in domain
        	    {
                	cellVol = mesh.V()[cellID];

                	//NP for 2 different radii
                	for(scalar r = 0.623926*radius;r < radius; r+=0.293976*radius )
                	{
                	    //NP try 8 subpoint derived from spherical coordinates
	                	for (scalar zeta=pi/4.;zeta<(2.*pi);zeta+=(pi/2.))
	                	{
                        	    for (scalar theta=(pi/4.);theta<pi;theta+=(pi/2.))
	                	    {
	                        	offset[0]=double(r)*Foam::sin(theta)*Foam::cos(zeta);
	                        	offset[1]=double(r)*Foam::sin(theta)*Foam::sin(zeta);
	                        	offset[2]=double(r)*Foam::cos(theta);
                            		#include "setWeightedSource.H"   // set source terms at position+offset
	                	    }
                    		}
	                	//NP try 2 more subpoints for each coordinate direction (6 total)
	                	for (int j=-1;j<=1;j+=2)
	                	{
	    	        	    offset[0]=double(r)*(double(j));
	                	    offset[1]=double(0.);offset[2]=double(0.);
                        	    #include "setWeightedSource.H"   //NP set source terms at position+offset

				    offset[1]=double(r)*(double(j));
	                	    offset[0]=double(0.);offset[2]=double(0.);
                        	    #include "setWeightedSource.H"   //NP set source terms at position+offset

	                	    offset[2]=double(r)*(double(j));
	                	    offset[0]=double(0.);offset[1]=double(0.);
                        	    #include "setWeightedSource.H"   //NP set source terms at position+offset
	                	}
                	}// end loop radiivoidfractions

	        	if(cellsSet>29 || cellsSet<0)
                	{
	        	    Info << "ERROR  cellsSet =" << cellsSet << endl;
        		}

                	//NP set source for particle center; source 1/nPts+weight of all subpoints that have not been found
                	scalar centreWeight = 1./nPoints*(nPoints-cellsSet);

                	// update voidfraction for each particle read
                	scalar newAlpha = mappedVoidfraction[cellID]- volume*centreWeight/cellVol;
                	if(newAlpha > alphaMin) mappedVoidfraction[cellID] = newAlpha;
                	else
                	{
                	    mappedVoidfraction[cellID] = alphaMin;
                	    tooMuch += (alphaMin-newAlpha) * cellVol;
                	}

                	// store cellweight for each particle --- this should be done for subpoints as well!!
                	particleWeights[index][0] += centreWeight;

                	// store particleVolume for each particle
                	particleVolumes[index][0] += volume*centreWeight;

        	    }// end if in cell
        	//NP reset counter of lost volume
        	if(index == nP[me]-1) Info << tab << "Total particle volume neglected: " << tooMuch << endl;
	    }// end loop all particles
    }
    else
    { 	
	    for(label index = 0; index < nP[me]; index++)
	    {

        	    cellsPerParticle[index] = 1;
        	    cellID = cellIDs[index];

        	    if (cellID > -1)  
        	    {
                	cellVol = mesh.V()[cellID];
        	    	radius = radii[me][index];
        	    	scalar volume =  4./3.*pi*radius*radius*radius;;

                	// Cell centered weight
                	particleWeights[index][0] = 1.;
			
			// store particleVolume for each particle
                	particleVolumes[index][0] = volume;

                	mappedVoidfraction[cellID] -= volume/cellVol;
        	    }
	    }
    }
    
    Info << tab << "Writing mapped volume fraction... " << endl;
    mappedVoidfraction.write();
    // Update mappedVoidfraction
    //forAll(mesh.cells(),index) mappedVoidfraction[index] = 1.;  

}
