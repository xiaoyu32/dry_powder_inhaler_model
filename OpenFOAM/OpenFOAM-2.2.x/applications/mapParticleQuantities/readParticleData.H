		 		  
//- We do not keep particle data @t=0
if ( runTime.timeName() != "0" )
{
	int count = runTime.value() / DEMts;				
	//Read only if master
	if(Pstream::master())
	{		
		Info<< tab << "Reading particle variables" << endl;
		// get path to particle VTK files
		string HH = string(filename);
		const char* charFilename;
		charFilename = HH.c_str();
		char index[100];
		sprintf(index, charFilename, count);

		fileName H(relativePath/index);
		Info << tab << "Opening file: " << H << endl;

		// set file pointer
		HH = string(H);
		const char * particleFilePath = HH.c_str();
		ifstream inputPtr(particleFilePath);

		// Read data
		string just_read = " ";
		int just_read_int;
		scalar just_read_scalar;
		
		while( getline( inputPtr, just_read ) )
		{		    
		    if( strcmp( just_read.c_str(), "ITEM: NUMBER OF ATOMS" ) == 0 )
		    {		        
			getline( inputPtr, just_read );			
			nPAll = atoi( just_read.c_str() );
			break;				
		    }		    		    
		}
		
		//nPAll = 10000000;
		Info << tab << "Number of particles = " << nPAll << endl;
				
		scalar fieldDummy;		
		while( getline( inputPtr, just_read ) )
		{		    		    
		    if( strncmp( just_read.c_str(), "ITEM: ATOMS id type", 19 ) == 0 )
		    {		       
			for(int index = 0; index < nPAll; ++index)
			{			     
			     
			     if ( index%1000000 ==0 ) Info << tab << "Just read particle " << index << "th." << endl;
			     inputPtr >> fieldDummy
				      >> fieldDummy
				      >> partAllRadius
				      >> partAllPos[0]
				      >> partAllPos[1]
				      >> partAllPos[2]
				      >> fieldDummy
				      >> partAllVel[0]
				      >> partAllVel[1]
				      >> partAllVel[2]
				      >> fieldDummy
				      >> fieldDummy
				      >> fieldDummy
				      >> fieldDummy
				      >> fieldDummy
				      >> fieldDummy;
			    			     	
			    //- Locate particles into CPU domains
			    if( Pstream::parRun() )
			    {
			      for(int ic = 0; ic < ncpus; ic++)
    			      {
				   int iCell = 0; 
				   for(int idir=0; idir<3; idir++)
				   {				 
				      scalar iCell_lo = ( partAllPos[idir] - coords[6*ic+2*idir] );
				      scalar iCell_hi = ( partAllPos[idir] - coords[6*ic+2*idir+1]);
				       if( iCell_lo >= 0 && iCell_hi <= 0 )
					   iCell++;
				   }
				   if ( iCell == 3 )
				   {
			    	     nP[ic]++;
				     radii[ic].append(partAllRadius);
				     positions [ic].append(partAllPos);
				     velocities[ic].append(partAllVel);
				   }  
			      }
			    }else
			    {
			       nP[me]++;
			       radii[me].append(partAllRadius);
			       positions [me].append(partAllPos);
			       velocities[me].append(partAllVel);
			    }   			    			    
			}			
			break;			
		    }	    
		}				
		
	       Info << tab << "Reading from the file is done..." << endl; 
	       
 	}
	
	
	if(Pstream::master())
	{
          for(int slave=Pstream::firstSlave(); slave<=Pstream::lastSlave(); slave++ )
          {
            Perr << tab << "master sending to slave " << slave << endl;
            OPstream toSlave(Pstream::scheduled, slave);
            toSlave << nP[slave];

	  }	  
	}else
	{
	     //- Receive data
            Perr<< tab << "slave receiving from master "
                    << Pstream::masterNo() << endl;
            IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());
            fromMaster >> nP[me];		
	}

        Pout << tab << "Local number of particles = "  << nP[me] << endl;
        Info << tab << "Number of particles per CPU was sent " << endl;

	if(Pstream::master())
	{	
            for(int slave=Pstream::firstSlave(); slave<=Pstream::lastSlave(); slave++ )
            {	
		OPstream toSlaveRad(Pstream::scheduled, slave);
        	toSlaveRad << radii[slave]; 
	    }
	}else
	{
	    IPstream fromMasterRad(Pstream::scheduled, Pstream::masterNo());
	    fromMasterRad >> radii[me];
	}
		
	Info << tab << "Radii were distibuted " << endl;

	if(Pstream::master())
	{	
            for(int slave=Pstream::firstSlave(); slave<=Pstream::lastSlave(); slave++ )
            {				
        	OPstream toSlavePos(Pstream::scheduled, slave);
        	toSlavePos << positions[slave]; 
            }
	}else  
	{
	    IPstream fromMasterPos(Pstream::scheduled, Pstream::masterNo());
	    fromMasterPos >> positions[me];		
	}	

	Info << tab << "Positions were distibuted " << endl;
	
	if(Pstream::master())
	{
            for(int slave=Pstream::firstSlave(); slave<=Pstream::lastSlave(); slave++ )
            {
        	OPstream toSlaveVel(Pstream::scheduled, slave);
        	toSlaveVel << velocities[slave]; 		   		   		   
            }
	}else  
	{
	    IPstream fromMasterVel(Pstream::scheduled, Pstream::masterNo());					
	    fromMasterVel >> velocities[me];					
	}
	
	Info << tab << "Velocities were distributed " << endl;
    			    
	//-Resize cellIDs
	cellIDs.setSize(nP[me],-1);	
	vector position(0,0,0);		
	for(int index = 0; index < nP[me]; ++index)
	{
	  //Pout << " index = " << index << " radii[me][index] " << radii[me][index] << endl;
	  if( radii[me][index] > SMALL)
	  {
              // create position vector
              for(int i=0;i<3;i++) position[i] = positions[me][index][i];
              // find cell
              cellIDs[index] = searchEngine.findCell(position,cellIDs[index],treeSearch);
	  }else{
	      Info<<"Small particles!"<<endl;
	  }
	}
				    
	if(verbose)
	{
		for(label ii=0; ii < exList.size(); ii++)
		{
			Info << "" << endl;
			int exIndex = exList[ii]; 
			{					
				Pout << tab << "index  = " << exIndex             << endl;
				Pout << tab << "rp     = " << radii[me][exIndex]      << endl;
				//Pout << tab << "Type   = " << types[exIndex]      << endl;
				Pout << tab << "Vp     = " << velocities[me][exIndex] << endl;
				//Pout << tab << "Omegap = " << omegas[exIndex]     << endl;			
				Pout << tab << "Xp     = " << positions[me][exIndex]  << endl;
				Pout << tab << "CellID = " << cellIDs[exIndex]    << endl;
				Pout << "" << endl;
			}	
		}
	}

}
			
