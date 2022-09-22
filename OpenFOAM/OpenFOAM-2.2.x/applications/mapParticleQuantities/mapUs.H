
{
    weightField.internalField() = 0.;
    label cellI(-1);
    vector valueVec(0,0,0);
    scalar weightP(0.);
    
    bool debug(false);

    for(int index=0; index < nP[me]; index++)
    {
	// Centering cell
	cellI = cellIDs[index];
        if (cellI > -1)
        {
            for(int iDir=0; iDir<3; iDir++) valueVec[iDir] = velocities[me][index][iDir];
            weightP = particleWeights[index][0];	    
            // first entry in this cell
            if(weightField[cellI] == 0)
            {
        	mappedUs[cellI] = valueVec;
        	weightField[cellI] = weightP;
            }
            else
            {
        	mappedUs[cellI] = (mappedUs[cellI]*weightField[cellI]+valueVec*weightP)/(weightField[cellI]+weightP);
        	weightField[cellI] += weightP;
            }
        }
	
	if(debug) Info << " index = " << index << " cellI = " << cellI << " weightP = " << weightP;
	
	// Neighboring cells    
	for(int subCell=1 ; subCell < maxCellPerPart; subCell++)
	{
            cellI = neighboringCellIDs[index][subCell];
            if (cellI > -1)
            {
        	for(int iDir=0; iDir<3; iDir++) valueVec[iDir] = velocities[me][index][iDir];
        	weightP = particleWeights[index][subCell];
        	// first entry in this cell
        	if(weightField[cellI] == 0)
        	{
        	    mappedUs[cellI] = valueVec;
        	    weightField[cellI] = weightP;
        	}
        	else
        	{
        	    mappedUs[cellI] = (mappedUs[cellI]*weightField[cellI]+valueVec*weightP)/(weightField[cellI]+weightP);
        	    weightField[cellI] += weightP;
        	}
 	        if(debug) Info << " weight2 = " << weightP;
            }
	}
	if(debug) Info << endl;	
    }
    Info << tab << "Writing mapped Us... " << endl;
    mappedUs.write();
}    
