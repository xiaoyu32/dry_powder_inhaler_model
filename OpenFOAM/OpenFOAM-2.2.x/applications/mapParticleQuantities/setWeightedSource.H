
// set force and mass source at cell subCellID
{
    // locate subPoint
    vector subPosition = position + offset;
    label partCellId(-1);
    partCellId = searchEngine.findCell(subPosition,partCellId,treeSearch); 
    
    if( partCellId != cellID )
    {
        nUnEqual++;
    }
    if( partCellId < 0 )
    {
        nNotFound++;
    }
    nTotal++;

    if (partCellId > -1)  // subPoint is in domain
    {
        // update voidfraction for each particle read
        scalar partCellVol = mesh.V()[partCellId];;// U.mesh().V()[partCellId];
        scalar particleVolume = volume/static_cast<scalar>(nPoints);
        scalar newAlpha = mappedVoidfraction[partCellId] - particleVolume / partCellVol;
        
        if( newAlpha > alphaMin ) mappedVoidfraction[partCellId] = newAlpha;
        else
        {
            mappedVoidfraction[partCellId] = alphaMin;
            tooMuch += (alphaMin-newAlpha) * partCellVol;
        }

        cellsSet++; // inc counter of cells affected

        //====================================================//
        // add sub particle representation
        bool createNew = true;
        label storeInIndex=0;
	
        for(int i=0; i < cellsPerParticle[index]; i++)
        {
            if( partCellId == cellIDs[index] ) 
            {
                storeInIndex = i;
                createNew = false;
                break;
            }else if( partCellId == neighboringCellIDs[index][i] )
	    {
                storeInIndex = i;
                createNew = false;
                break;	    
	    }
        }	    	

        if(createNew)
        {
            cellsPerParticle[index] ++;
            storeInIndex = cellsPerParticle[index]-1;
            neighboringCellIDs[index][storeInIndex] = partCellId;
        }

        particleWeights[index][storeInIndex] += 1/static_cast<scalar>(nPoints);
        particleVolumes[index][storeInIndex] += particleVolume;
        //====================================================//
    }
}
