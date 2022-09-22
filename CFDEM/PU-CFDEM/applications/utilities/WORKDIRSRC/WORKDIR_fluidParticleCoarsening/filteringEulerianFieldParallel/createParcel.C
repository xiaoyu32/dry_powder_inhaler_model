#include "createParcel.H"
 
namespace Foam
{ 
 	void createParcels
	(
		const fvMesh& mesh,			
		cfdemCloud& sm,
		cfdemCloud& parcelsm,
		const int& nparticle_,
		const bool verbose_,
		const fileName outputRelativePath_, // Write parcel diameter in a file in the directory of ./PostProcessing..	
		labelListList& parcelList,
		const bool& weighting_,
		const scalar& minX_,
		const scalar& maxX_,
		const scalar& minY_,
		const scalar& maxY_,
		const scalar& minZ_,
		const scalar& maxZ_,
		const labelList& partIDInSubVolume_,
		labelList& parcelIDInSubVolume_					 			
	)
	{		

		//if ( nparticle_ > sm.numberOfParticles() )
		if ( nparticle_ > partIDInSubVolume_.size() )
		{
			FatalError << " Number of particles in a parcel > number of particles" << abort(FatalError);
		}
		
		if ( nparticle_ < 1 )
		{
			FatalError << " Number of particles < 0 in a parcel " << abort(FatalError);
		}
		
		// Number of particles in a parcel
		int k = nparticle_; 		
				
		// Dimensions, exact OR approximate  
		int dim =3; double eps = 0;
						
		// Number of points
		int nPts;
		//nPts =  sm.numberOfParticles();		
		nPts = partIDInSubVolume_.size();
		
		Pout << " Number of particles in a parcel = " << nparticle_ << endl;		
		
		// Squared radius
		const scalar sqRad =  ( (maxX_ - minX_) * (maxX_ - minX_) )
				     +( (maxY_ - minY_) * (maxY_ - minY_) )
				     +( (maxZ_ - minZ_) * (maxZ_ - minZ_) );			
        /*              
        const scalar sqRad = ( ( (maxX_ - min_x) * (maxX_ - min_x) )
                              +( (max_y - min_y) * (max_y - min_y) )
                              +( (max_z - min_z) * (max_z - min_z) ) ) / 10.;  
		*/  
		
		Pout << " Squared radius = " << sqRad << endl;
		
		// Data points
		ANNpointArray dataPts;
		// Query points
		ANNpoint queryPt;
		
		ANNidxArray	nnIdx;			// 	near neighbour indices
		ANNdistArray 	dists;			//	near neighbour distances
		ANNkd_tree* 	kdTree;			//	search structure
		
		// Allocate 
		queryPt = annAllocPt(dim);
		dataPts = annAllocPts(nPts, dim);
		nnIdx = new ANNidx[k];
		dists = new ANNdist[k];		
		
		// Create particle list
		//labelList particleList(sm.numberOfParticles());				
		labelList particleList(partIDInSubVolume_.size());
		
		//for(int index = 0; index <  sm.numberOfParticles(); index++)
	    //    {			
	    for(int ii =0; ii < partIDInSubVolume_.size(); ii++)
		{
				int index = partIDInSubVolume_[ii];
				//dataPts[index][0] = sm.position(index).x();		
				//dataPts[index][1] = sm.position(index).y();		
				//dataPts[index][2] = sm.position(index).z();
				//particleList[index] = index;			
				dataPts[ii][0] = sm.position(index).x();		
				dataPts[ii][1] = sm.position(index).y();		
				dataPts[ii][2] = sm.position(index).z();
				particleList[ii] = index;					
		}	
		
						
		kdTree = new ANNkd_tree(dataPts, nPts, dim);
				
		label parcelI = 0;						
		forAll(particleList,ii)		
		{				
						
			if ( particleList[ii] > -1 )
				{					
					
					label partGlobalID = particleList[ii];
					//if (verbose_) Pout << " Parcel = " << parcelI << endl;
					
					//queryPt[0] = sm.position(index).x();
					//queryPt[1] = sm.position(index).y();
					//queryPt[2] = sm.position(index).z();
										
					queryPt[0] = sm.position(partGlobalID).x();
					queryPt[1] = sm.position(partGlobalID).y();
					queryPt[2] = sm.position(partGlobalID).z();
	
					kdTree->annkFRSearch(
					                		queryPt,			// query point					
											sqRad,				// squared radius
											k,				// number of the near neighbours to return
											nnIdx,				// nearest neighbor array
											dists,				// dist to near neighbours
											eps			);					
							
					//if (verbose_) Pout << "NN:  Index  Distance\n";
																				
					int sumnnIdx = 0; 
					 	
					for (int i = 0; i < k; i++)
					{						
						if ( particleList[nnIdx[i]] == -1 )
						{
							
							  //particleList[nnIdx[i]] = sm.numberOfParticles() + 1;
							particleList[nnIdx[i]] = partIDInSubVolume_.size() + 1;
								       sumnnIdx += particleList[nnIdx[i]];
									   
						}
						else
						{
							//parcelList[parcelI][i] = nnIdx[i] ;
							  parcelList[parcelI][i] = particleList[nnIdx[i]];
							  
						   	  //Pout << " nnIdx " << nnIdx[i] << " particleList " << particleList[nnIdx[i]] << endl;
							  
							  particleList[nnIdx[i]] = -1 ;
								           sumnnIdx += particleList[nnIdx[i]];
										   
							
						}		 						 
					}
										
					if ( sumnnIdx == -nparticle_ )
					{
						parcelI++;																				
					}
					
				}		
		}	
								
		// Local number of parcel 
		int npLocalParcel = parcelI;
								
		// Parallel computation
		int totalNpParcel = parcelI;
		if(Pstream::parRun())
		{
			reduce(totalNpParcel, sumOp<label>());
		}	
		
		Pout << " Number of parcels = " << npLocalParcel << endl;
		//Pout << " Number of particles not used to construct parcels = " << sm.numberOfParticles()-parcelI*nparticle_ << endl;
		Pout << " Number of particles not used to construct parcels = " << partIDInSubVolume_.size() - npLocalParcel*nparticle_ << endl;
						
		// Set number of parcels
		//parcelsm.setNumberOfParticles(totalNpParcel);	
		parcelsm.setNumberOfParticles(totalNpParcel);	
								
		// Resize parcel list 
		//parcelList.resize(parcelI);	
		parcelList.resize(npLocalParcel);
		
		// Parcel center of mass
		//vectorField parcelCenter(parcelI,vector(0,0,0));
		vectorField parcelCenter(npLocalParcel,vector(0,0,0));
		// Calculate weighted velocities
		//vectorField parcelVel(parcelI,vector(0,0,0));
		vectorField parcelVel(npLocalParcel,vector(0,0,0));
		// Parcel diameter
		//scalarField parcelDiameter(parcelI,scalar(0.0));
		scalarField parcelDiameter(npLocalParcel,scalar(0.0));	
								
		vector Up(0,0,0);	
		vector position(0,0,0);						
		//for(int parcelII = 0; parcelII < parcelsm.numberOfParticles(); parcelII++)		
		for(int parcelII = 0; parcelII < npLocalParcel ; parcelII++)
		{ 			
			labelList parcelListL = parcelList[parcelII];
			Up = vector(0,0,0);
			position = vector(0,0,0);
			parcelDiameter[parcelII] = 0;
						
			forAll(parcelListL,particleI)
			{				
					
				Up = sm.velocity(parcelListL[particleI]);
				position = sm.position(parcelListL[particleI]);	
				 					
				for(int j=0;j<3;j++) 
				{	
						// Parcel velocity 
						parcelVel[parcelII][j]    += Up[j];
						// Parcel center of mass
						parcelCenter[parcelII][j] += position[j];
		        	      					
				}	
				parcelDiameter[parcelII] =  2.*nparticle_*sm.radius(parcelListL[particleI]);											
			}				
		}	
								
		List < scalar > listScalar(1);	 
		List < List < scalar > > listScalarScalar(npLocalParcel,listScalar);
		List < List < List < scalar > > > listParcelRadiis(Pstream::nProcs(),listScalarScalar);
		
		List < vector > listVector(1);	 
		List < List < vector > > listVectorVector(npLocalParcel,listVector);
		List < List < List < vector > > > listParcelPositions(Pstream::nProcs(),listVectorVector);					
		List < List < List < vector > > > listParcelVelocities(Pstream::nProcs(),listVectorVector);	
			
		for(int index = 0; index <  npLocalParcel; index++)
		{
			
			listParcelRadiis[Pstream::myProcNo()][index][0] = parcelDiameter[index]/2.;

			listParcelPositions[Pstream::myProcNo()][index] = vector(parcelCenter[index][0]/nparticle_,
														 		     parcelCenter[index][1]/nparticle_,
														 			 parcelCenter[index][2]/nparticle_);

 			listParcelVelocities[Pstream::myProcNo()][index] = vector(parcelVel[index][0],
 														 		      parcelVel[index][1],
 																      parcelVel[index][2]);	
																	  
		}
		
		double **parcelPositions;
		double **parcelVelocities;
		double **parcelRadii;
		double **parcelCellID;	

		parcelsm.dataExchangeM().allocateArray(parcelPositions,0.,3);
		parcelsm.dataExchangeM().allocateArray(parcelVelocities,0.,3);
		parcelsm.get_radii(parcelRadii); 
		parcelsm.get_cellIDs(parcelCellID);	
				

		// Distrubute to all processors
		Pstream::gatherList(listParcelRadiis);
		Pstream::scatterList(listParcelRadiis);

		Pstream::gatherList(listParcelPositions);
		Pstream::scatterList(listParcelPositions);

		Pstream::gatherList(listParcelVelocities);
		Pstream::scatterList(listParcelVelocities);			
											
		// Create global parcels
		List < List < scalar > > globalListParcelRadiis; //(totalNpParcel); //(1,-1);
	
		globalListParcelRadiis = ListListOps::combine< List < List < scalar > > >
		(
			listParcelRadiis,
			accessOp< List < List < scalar > > >()
		);			
		
		// Create global parcels
		List < List < vector > > globalListParcelPositions; //(totalNpParcel); //(1,-1);

		globalListParcelPositions = ListListOps::combine< List < List < vector > > >
		(
			listParcelPositions,
			accessOp< List < List < vector > > >()
		);			
					
		// Create global parcels
		List < List < vector > > globalListParcelVelocities; //(totalNpParcel); //(1,-1);

		globalListParcelVelocities = ListListOps::combine< List < List < vector > > >
		(
			listParcelVelocities,
			accessOp< List < List < vector > > >()
		);			
		
		// Define Parcel Cloud positions, velocities, radius and cellIDs				
	
		if(Pstream::master)
		{	

			for(int index = 0; index <  parcelsm.numberOfParticles(); index++)
			{	
				for(int j=0;j<3;j++) 
				{			
					parcelPositions[index][j] = globalListParcelPositions[index][0][j];
					parcelVelocities[index][j] = globalListParcelVelocities[index][0][j] ;		
				}							
				parcelRadii[index][0] = globalListParcelRadiis[index][0] ;			

			}
			
		}
				
		// Send to child CPUs
		if(Pstream::parRun())
		{
			//for(int index = 0; index < maxNumberOfParticles; index++)
			for(int index = 0; index < parcelsm.numberOfParticles(); index++)
			{
				for(int idir = 0; idir < 3; idir++)
				{
					Pstream::scatter(parcelPositions[index][idir]);
					Pstream::scatter(parcelVelocities[index][idir]);
				}
				Pstream::scatter(parcelRadii[index][0]);
			}
		}
							
		parcelsm.locateM().findCell(NULL,parcelPositions,parcelCellID,parcelsm.numberOfParticles());
		parcelsm.setPos(parcelPositions);
		parcelsm.setVel(parcelVelocities);
		
		// Parcel IDs in sub-volume
		int ii = 0 ;
		position = vector(0,0,0);
		
		for(int index = 0; index < parcelsm.numberOfParticles(); index++)
		{
			label cellI = parcelsm.cellIDs()[index][0];						
			if(cellI > -1)
			{						      
				parcelIDInSubVolume_[ii] = index; 
				ii++;	            							      
			}
		}
		parcelIDInSubVolume_.resize(ii);
					
		// Number of parcell WARNING There is a problem here Why parcelIDInSubVolume_.size() != npLocalParcel
		int trimNpParcel = pow(10,10);		
		if( parcelIDInSubVolume_.size() < trimNpParcel ) trimNpParcel = parcelIDInSubVolume_.size();		
		if( npLocalParcel < trimNpParcel ) trimNpParcel = npLocalParcel;
															
		Pout << " Create parcel-> Trimming of max. number of parcel: " << trimNpParcel 
			 << " parcel sub-volume list: " << parcelIDInSubVolume_.size() << " parcel list size: " << npLocalParcel << endl;	
		/*												
		if( verbose_ )
		{	
			//for(int parcelII = 0; parcelII < parcelsm.numberOfParticles(); parcelII++)
			//for(int parcelII = 0; parcelII < parcelIDInSubVolume_.size(); parcelII++)
			for(int parcelII = 0; parcelII < trimNpParcel; parcelII++)
			{ 
				Pout << " ParcelID = "  << parcelII << " Cell ID " << parcelsm.cellIDs()[parcelIDInSubVolume_[parcelII]][0] 
					                    << " Position " << parcelsm.position(parcelIDInSubVolume_[parcelII])
									    << " Particle positions ";
										for(int ii = 0; ii < parcelList[parcelII].size(); ii++ )
										{		
											Pout << sm.position(parcelList[parcelII][ii]);
										}		 
										Pout << endl;
			}
		}
		*/			
				
		if ( verbose_ && nparticle_ <=2)
		{							
			int index = 0;
			labelList parcelPart = parcelList[index];
			Pout << " Parcel particle list " << parcelPart << endl; 
			Pout << " Parcel center     " <<  parcelsm.position(parcelIDInSubVolume_[index])  		<< endl;	
			Pout << " Parcel velocity   " <<  parcelsm.velocity(parcelIDInSubVolume_[index])  		<< endl;
			Pout << " Parcel diameter   " << 2.*parcelsm.radius(parcelIDInSubVolume_[index]) 		<< endl;
					
			forAll(parcelPart, ii)
			{
				Pout << " Particle " << parcelPart[ii] << endl;
				Pout << " Particle center    " <<        sm.position(parcelPart[ii])   << endl;
				Pout << " Particle velocity  " <<        sm.velocity(parcelPart[ii])   << endl;
				Pout << " Particle diameter  " <<       2.*sm.radius(parcelPart[ii]) 	<< endl;
			}			
		}
												
	}
 

}	
