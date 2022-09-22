#include "output.H"

namespace Foam
{

	// Initiate output file
	void initOutput
	(	
		const int&          nParticleClass_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		OFstream*&          outputCollFile_,
		const bool&             bidisperse_,
		const bool&          calcCollision_		
	)
	{
	
	    	outputKinFile_ =  new OFstream(mesh_.time().path()/outputRelativePath_/"kineticStress");
		*outputKinFile_ << "#Time  " << tab << " "
		                << "Vol    " <<        " "  
	                	<< "alpp   " << tab << " "
				<< "Kin_XX " << tab << " "
				<< "Kin_YY " << tab << " "
				<< "Kin_ZZ " << tab << " " 
				<< "Kin_XY " << tab << " " 
				<< "Kin_XZ " << tab << " " 
				<< "Kin_YZ " << tab << " "; 		

		if( bidisperse_ ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        	{
				*outputKinFile_  << "Kin_XX["<<iPartClass<<"] " << tab << " "
						<< "Kin_YY["<<iPartClass<<"] " << tab << " "
						<< "Kin_ZZ["<<iPartClass<<"] " << tab << " " 
						<< "Kin_XY["<<iPartClass<<"] " << tab << " " 
						<< "Kin_XZ["<<iPartClass<<"] " << tab << " " 
						<< "Kin_YZ["<<iPartClass<<"] " << tab << " "; 
			}
		}
		*outputKinFile_ << endl;
	
		if( calcCollision_ )
		{
    			outputCollFile_ =  new OFstream(mesh_.time().path()/outputRelativePath_/"collisonalStress");
			*outputCollFile_ << "#Time  " << tab << " "
		                	 << "Vol    " << 	" "  
	                		 << "alpp   " << tab << " "			 
	                		 << "Coll_XX " << tab << " "
					 << "Coll_YY " << tab << " "
					 << "Coll_ZZ " << tab << " " 
					 << "Coll_XY " << tab << " " 
					 << "Coll_XZ " << tab << " " 
					 << "Coll_YZ " << tab << " "; 	
			if( bidisperse_ ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        		{
					*outputCollFile_  << "Coll_XX["<<iPartClass<<"] " << tab << " "
							 << "Coll_YY["<<iPartClass<<"] " << tab << " "
							 << "Coll_ZZ["<<iPartClass<<"] " << tab << " " 
							 << "Coll_XY["<<iPartClass<<"] " << tab << " " 
							 << "Coll_XZ["<<iPartClass<<"] " << tab << " " 
							 << "Coll_YZ["<<iPartClass<<"] " << tab << " "; 			
				}
			}
			*outputCollFile_ << endl;	
		}
	}
	
	// Write into output file
	void writeOutput
	(	
		const int&          nParticleClass_,
		const Time& 		   runTime_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		OFstream*&          outputCollFile_,
		const bool&             bidisperse_,
		const bool&          calcCollision_,
		SymmTensor<double>*       sigmaKin_,					
		SymmTensor<double>*	 sigmaColl_,
		const scalar&		 domainVol_,
		const scalar&	   alppInSubVolume_	)
	{
	
		// Write output
		Info << " " << endl;
		Pout << " Writing particle kinetic stresses into the file " << mesh_.time().path()/outputRelativePath_/"kineticStress" << endl;

		int iPartClass = 0;
		*outputKinFile_ << runTime_.value() 		       << tab << " " 
		                << domainVol_ 			       << tab << " "
				<< alppInSubVolume_		       << tab << " "   	
				<< sigmaKin_[iPartClass][0]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][1]/domainVol_ << tab << " " 
				<< sigmaKin_[iPartClass][2]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][3]/domainVol_ << tab << " " 
				<< sigmaKin_[iPartClass][4]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][5]/domainVol_ << tab << " " ;								

		if( bidisperse_ ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        	{
				*outputKinFile_ << sigmaKin_[iPartClass][0]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][1]/domainVol_ << tab << " " 
						<< sigmaKin_[iPartClass][2]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][3]/domainVol_ << tab << " " 
						<< sigmaKin_[iPartClass][4]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][5]/domainVol_ << tab << " " ; 
			}
		}
		*outputKinFile_ << endl;

		if( calcCollision_ )
		{
    			Pout<< " Writing particle kinetic stresses into the file " << mesh_.time().path()/outputRelativePath_/"collisionalStress" << endl;
			iPartClass = 0;
			*outputCollFile_ << runTime_.value() 		         << tab << " " 	
					 << domainVol_ 			         << tab << " "
					 << alppInSubVolume_		         << tab << " "  
					 << sigmaColl_[iPartClass][0]/domainVol_ << tab << " " 		
					 << sigmaColl_[iPartClass][1]/domainVol_ << tab << " " 
					 << sigmaColl_[iPartClass][2]/domainVol_ << tab << " " 		
					 << sigmaColl_[iPartClass][3]/domainVol_ << tab << " " 
					 << sigmaColl_[iPartClass][4]/domainVol_ << tab << " " 		
					 << sigmaColl_[iPartClass][5]/domainVol_ << tab << " " ;
			if( bidisperse_ ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        		{
					*outputCollFile_ << sigmaColl_[iPartClass][0]/domainVol_ << tab << " " 		
					  		 << sigmaColl_[iPartClass][1]/domainVol_ << tab << " " 
							 << sigmaColl_[iPartClass][2]/domainVol_ << tab << " " 		
							 << sigmaColl_[iPartClass][3]/domainVol_ << tab << " " 
							 << sigmaColl_[iPartClass][4]/domainVol_ << tab << " " 		
							 << sigmaColl_[iPartClass][5]/domainVol_ << tab << " " ;			
				}
			}
			*outputCollFile_ << endl;	
		}	
	}	
	
}	
