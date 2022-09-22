//write particle information (works only in serial)

{

	List<scalar> partVolFractions(nPAll);	
	List<scalar> contactNumber(nPAll);

	for( label index = 0; index < nPAll; ++index )
	    contactNumber[index] = 0.0; 

	for( label index = 0; index < nPAll; ++index )
	{
	    
	    label cellID = cellIDs[index];
	    
	    if( cellID > 0 )
	       partVolFractions[index] = mappedVoidfraction[cellID];
	    else
	    {
	       partVolFractions[index] = 0.0;
	       Info<<"Warning no cell information!"<<index<<" "<<cellID<<endl;
	    }
	    
	    //compute the particle contact number
	    /*for( label ii = index + 1; ii < nPAll; ++ii )
	    {
	        
		//this is actually diameter
		if( Foam::mag( positions[0][index] - positions[0][ii] ) < partAllRadius )
		{
		    contactNumber[index] += 1.0;
		    contactNumber[ii] += 1.0;
		}
		
	    }*/
	    
	       
	}
	
	//write info
	ofstream outputStream( outputFilename );	
	
	//header
	outputStream<<"ITEM: TIMESTEP"<<endl;	
	outputStream<<"1000"<<endl;     			//hardcoding !!!
	outputStream<<"ITEM: NUMBER OF ATOMS"<<endl;
	outputStream<<nPAll<<endl;
	outputStream<<"ITEM: BOX BOUNDS ff pp pp"<<endl;
	outputStream<<"0 0.0045"<<endl; 			//hardcoding !!!
	outputStream<<"0 0.0045"<<endl; 			//hardcoding !!!
	outputStream<<"0 0.018"<<endl;  			//hardcoding !!!
	outputStream<<"ITEM: ATOMS id type diameter x y z vx vy vz volFrac nCont"<<endl;
	
	for( label index = 0; index < nPAll; ++index )
	{
	    
	    
	    outputStream<<index<<" "
	    	    	<<"1"<<" "
		   	<<partAllRadius<<" "
		    	<<positions[0][index].x()<<" "
		    	<<positions[0][index].y()<<" "
		    	<<positions[0][index].z()<<" "
		    	<<velocities[0][index].x()<<" "
		    	<<velocities[0][index].y()<<" "
		    	<<velocities[0][index].z()<<" "
		    	<<partVolFractions[index]<<" "
		    	<<contactNumber[index]<<endl;
		    
		    
	}	
	
	outputStream.close();
	
}
