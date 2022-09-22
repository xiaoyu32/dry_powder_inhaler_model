#include "createParcels.H"

namespace Foam
{
 	void createParcels
	(
		const fvMesh& mesh,			
		cfdemCloud& sm,
		const int& parcelSize_,		
		int**& parcelCloud_,
		double**& parcelPositions_,
		double**& parcelVelocities_,
		int*& parcelNparts_,				
		double** & parcelKinStress_,
		scalar& aveSubQparcel2_,
		vector& meanParcelVel_,
		const bool verbose_ 			
	)
	{		

		if ( parcelSize_ * parcelSize_ * parcelSize_  > sm.numberOfParticles() )
		{
			FatalError << " Number of particles in a parcel > number of particles" << abort(FatalError);
		}
		
		if ( parcelSize_ < 1 )
		{
			FatalError << " Number of particles < 0 in a parcel " << abort(FatalError);
		}
		
		// Number of particles in a parcel
		int k = parcelSize_ * parcelSize_ * parcelSize_; 		
				
		// Dimensions, exact OR approximate  
		int dim =3; double eps = 0;
						
		// Number of points
		int nPts;
		nPts =  sm.numberOfParticles();		
				
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
		 		
		for(int index = 0; index < sm.numberOfParticles(); index++)
	        {			
				dataPts[index][0] = sm.position(index).x();		
				dataPts[index][1] = sm.position(index).y();		
				dataPts[index][2] = sm.position(index).z();
		}
				
		kdTree = new ANNkd_tree(dataPts, nPts, dim);
		
		// Initialize sub-parcel agitation	
		aveSubQparcel2_ = 0.;
		
		// Initialize parcel velocity
		meanParcelVel_ = vector(0,0,0);	
		
		for(int index = 0; index <  sm.numberOfParticles(); index++)
	        {				
			
			// Particle neighbouring search distance
			scalar sqRad = parcelSize_ * sm.radius(index);		
									
			queryPt[0] = sm.position(index).x();
			queryPt[1] = sm.position(index).y();
			queryPt[2] = sm.position(index).z();

			kdTree->annkFRSearch(
					        queryPt,			// query point					
						sqRad,				// squared radius
						k,				// number of the near neighbours to return
						nnIdx,				// nearest neighbor array
						dists,				// dist to near neighbours
						eps			);					
 				
			int nParts = 0;
			scalar dist = 0;	
			
			// Initialize parcel velocities & positions & kinetic stresses
			for(int j=0;j<3;j++) 
			{
				parcelVelocities_[index][j] = 0.;
				 parcelPositions_[index][j] = 0.;
				 parcelKinStress_[index][j] = 0.;
			       parcelKinStress_[index][2*j] = 0.;
			}
			
			for (int i = 0; i < k; i++)
			{
				parcelCloud_[index][i] = nnIdx[i];
				
				dist = mag( sm.position(nnIdx[i]) - sm.position(index) ) - sqRad;
				if ( dist < SMALL ) 
				{				
					for(int j=0;j<3;j++) 
					{	
							// Parcel velocity 
							parcelVelocities_[index][j] += sm.velocity(nnIdx[i])[j];	
							// Parcel center of mass
							 parcelPositions_[index][j] += sm.position(nnIdx[i])[j];	

					}				
					nParts++;
				}							
			}

			for(int j=0;j<3;j++) parcelPositions_[index][j] /= nParts;			
			parcelNparts_[index] = nParts;	
																													
			// Parcel kinetic stresses
			for(int i = 0; i < parcelNparts_[index]; i++)
			{
				int particleID = parcelCloud_[index][i]; 
				
				// U'xU'x
				parcelKinStress_[index][0] +=  ( sm.velocity(particleID)[0] - parcelVelocities_[index][0] )   
						   	      *( sm.velocity(particleID)[0] - parcelVelocities_[index][0] );

				// U'yU'y
                                parcelKinStress_[index][1] +=  ( sm.velocity(particleID)[1] - parcelVelocities_[index][1] )
                                                              *( sm.velocity(particleID)[1] - parcelVelocities_[index][1] );

				// U'zU'z
                                parcelKinStress_[index][2] +=  ( sm.velocity(particleID)[2] - parcelVelocities_[index][2] )
                                                              *( sm.velocity(particleID)[2] - parcelVelocities_[index][2] );
 
				// U'xU'y
                                parcelKinStress_[index][3] +=  ( sm.velocity(particleID)[0] - parcelVelocities_[index][0] )
                                                              *( sm.velocity(particleID)[1] - parcelVelocities_[index][1] );

				// U'xU'z
                                parcelKinStress_[index][4] +=  ( sm.velocity(particleID)[0] - parcelVelocities_[index][0] )
                                                              *( sm.velocity(particleID)[2] - parcelVelocities_[index][2] );

				// U'yU'z
                                parcelKinStress_[index][5] +=  ( sm.velocity(particleID)[1] - parcelVelocities_[index][1] )
                                                              *( sm.velocity(particleID)[2] - parcelVelocities_[index][2] );
			}
			
			// Mean parcel velocity
			for(int j=0;j<3;j++) meanParcelVel_[j] += parcelVelocities_[index][j];

			// Domain-averaged parcel agitation
			aveSubQparcel2_ += 1./2. * ( parcelKinStress_[index][0] + parcelKinStress_[index][1] + parcelKinStress_[index][2] );

		}
		
		for(int j=0;j<3;j++) meanParcelVel_[j] /= sm.numberOfParticles();
		
		if ( verbose_ )
		{							
			int index = 0;
			Info << " Parcel particle list ";
			for (int i = 0; i < parcelNparts_[index]; i++) 
			{
				Info << parcelCloud_[index][i] << " " ;
			}
			Info << endl;
			
			Info << " Parcel center     " <<  parcelPositions_[index][0] << "," <<  parcelPositions_[index][1] << "," <<  parcelPositions_[index][2] << endl;	
			Info << " Parcel velocity   " << parcelVelocities_[index][0] << "," << parcelVelocities_[index][1] << "," << parcelVelocities_[index][2] << endl;
					
			for (int i = 0; i < parcelNparts_[index]; i++)			
			{
				Info << " Particle " << parcelCloud_[index][i] << endl;
				Info << " Particle center    " <<        sm.position(parcelCloud_[index][i])   << endl;
				Info << " Particle velocity  " <<        sm.velocity(parcelCloud_[index][i])   << endl;
			}			
		}
																
	}

}
