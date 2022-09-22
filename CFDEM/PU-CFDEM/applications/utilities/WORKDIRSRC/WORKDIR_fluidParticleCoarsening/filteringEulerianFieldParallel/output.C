#include "output.H"

namespace Foam
{
	void writeBins
	(
		const Time& runTime,
	    	const fvMesh& mesh,
		const int& nBin_,
		const scalar& minalpp_,
		const scalar& maxalpp_,
		const scalar& minVelSlip_,
		const scalar& maxVelSlip_,
		const scalar& minTaupVelSlip_,
		const scalar& maxTaupVelSlip_,
		const fileName outputRelativePath_,		
		const int& minfilterwidth,
		const int& maxfilterwidth,
		const int& filterincrement,
		const int& nparticle_,
		Field< Field < Field <scalar> > >& 	giiCondalppFilter_,
		Field< Field < Field <scalar> > >& 	velSlipCondalppFilter_,
		Field < Field <scalar> > & 		taupVelSlipCondalppFilter_,
		Field < Field <scalar> > & 		numberCondalppFilter_,	
		Field< Field < Field <scalar> > >&	VelSlipJointgiiCondalppFilter_,
		Field< Field < Field <scalar> > >&	TaupVelSlipJointgiiCondalppFilter_,	
		Field <Field < Field <scalar> > >& 	numberVelSlipJointgiiCondalppFilter_,
		Field <Field < Field <scalar> > >& 	numberTaupVelSlipJointgiiCondalppFilter_,				
		bool verbose_,
		bool EulerEulerFiltering_,
		bool EulerianVslipBin_ 
	)		
	{

		// Create output folder	
		if( !isDir(mesh.time().path()/outputRelativePath_/runTime.timeName()) )
		{
			mkDir(mesh.time().path()/outputRelativePath_/runTime.timeName() );
		}
				
		// Filename
		char charfilterwidth[100]; 
		char charfilterwidthVelSlip[100]; 
		char charfilterwidthTaupVelSlip[100]; 		
		char charfilterwidthN[100]; 

		char charfilterwidthJointVelSlip[100]; 
		char charfilterwidthJointTaupVelSlip[100]; 		
		char charfilterwidthJointNVelSlip[100]; 
		char charfilterwidthJointNTaupVelSlip[100];
				
		// Debugging for this subroutine
		// verbose_ = true;
								
		for( int filterwidth = minfilterwidth; filterwidth <= maxfilterwidth; filterwidth += filterincrement )
		{

			int Filter = 2*(filterwidth+1)+1;

			sprintf(charfilterwidth,        	 "h_vs_baralpp_%dX%dX%d"			,Filter,Filter,Filter);			
			sprintf(charfilterwidthVelSlip, 	 "velSlip_vs_baralpp_%dX%dX%d"			,Filter,Filter,Filter);
			sprintf(charfilterwidthTaupVelSlip, 	 "TaupVelSlip_vs_baralpp_%dX%dX%d"		,Filter,Filter,Filter);
			sprintf(charfilterwidthN, 		 "Nrealization_%dX%dX%d"			,Filter,Filter,Filter);

			sprintf(charfilterwidthJointVelSlip, 	 "JointH_VelSlip_vs_baralpp_%dX%dX%d"		,Filter,Filter,Filter);
			sprintf(charfilterwidthJointTaupVelSlip, "JointH_TaupVelSlip_vs_baralpp_%dX%dX%d"	,Filter,Filter,Filter);
			sprintf(charfilterwidthJointNVelSlip,	 "JointNreal_VelSlip_%dX%dX%d"			,Filter,Filter,Filter);
			sprintf(charfilterwidthJointNTaupVelSlip,"JointNreal_TaupVelSlip_%dX%dX%d"		,Filter,Filter,Filter);

			
			fileName outputfile(charfilterwidth);
			fileName outputfileVelSlip(charfilterwidthVelSlip);			
			fileName outputfileTaupVelSlip(charfilterwidthTaupVelSlip);				
			fileName outputfileN(charfilterwidthN);
			
			fileName outputfileJointVelSlip(charfilterwidthJointVelSlip);			
			fileName outputfileJointTaupVelSlip(charfilterwidthJointTaupVelSlip);				
			fileName outputfileJointNVelSlip(charfilterwidthJointNVelSlip);			
			fileName outputfileJointNTaupVelSlip(charfilterwidthJointNTaupVelSlip);							
			
			if ( nparticle_ > 0 && !EulerianVslipBin_)
			{								
				sprintf(charfilterwidth,		 "h_vs_baralpp_%dX%dX%d_npart=%d"			,Filter,Filter,Filter,nparticle_);
				sprintf(charfilterwidthVelSlip,		 "velSlip_vs_baralpp_%dX%dX%d_npart=%d"			,Filter,Filter,Filter,nparticle_);
				sprintf(charfilterwidthTaupVelSlip,	 "TaupVelSlip_vs_baralpp_%dX%dX%d_npart=%d"		,Filter,Filter,Filter,nparticle_);
				sprintf(charfilterwidthN,		 "Nrealization_%dX%dX%d_npart=%d"			,Filter,Filter,Filter,nparticle_);

				sprintf(charfilterwidthJointVelSlip, 	 "JointH_VelSlip_vs_baralpp_%dX%dX%d_npart=%d"		,Filter,Filter,Filter,nparticle_);
				sprintf(charfilterwidthJointTaupVelSlip, "JointH_TaupVelSlip_vs_baralpp_%dX%dX%d_npart=%d"	,Filter,Filter,Filter,nparticle_);
				sprintf(charfilterwidthJointNVelSlip,	 "JointNreal_VelSlip_%dX%dX%d_npart=%d"			,Filter,Filter,Filter,nparticle_);
				sprintf(charfilterwidthJointNTaupVelSlip,"JointNreal_TaupVelSlip_%dX%dX%d_npart=%d"		,Filter,Filter,Filter,nparticle_);
				
				outputfile = charfilterwidth;
				outputfileVelSlip = charfilterwidthVelSlip;	
				outputfileTaupVelSlip = charfilterwidthTaupVelSlip;				
				outputfileN = charfilterwidthN;		
				
				outputfileJointVelSlip = charfilterwidthJointVelSlip;			
				outputfileJointTaupVelSlip = charfilterwidthJointTaupVelSlip;				
				outputfileJointNVelSlip = charfilterwidthJointNVelSlip;			
				outputfileJointNTaupVelSlip = charfilterwidthJointNTaupVelSlip;					
						
			}
			else if (EulerEulerFiltering_)
			{
				sprintf(charfilterwidth,        	 "Eulerian_h_vs_baralpp_%dX%dX%d"			,Filter,Filter,Filter);			
				sprintf(charfilterwidthVelSlip, 	 "Eulerian_velSlip_vs_baralpp_%dX%dX%d"			,Filter,Filter,Filter);
				sprintf(charfilterwidthTaupVelSlip, 	 "Eulerian_TaupVelSlip_vs_baralpp_%dX%dX%d"		,Filter,Filter,Filter);
				sprintf(charfilterwidthN, 		 "Eulerian_Nrealization_%dX%dX%d"			,Filter,Filter,Filter);

				sprintf(charfilterwidthJointVelSlip, 	 "Eulerian_JointH_VelSlip_vs_baralpp_%dX%dX%d"		,Filter,Filter,Filter);
				sprintf(charfilterwidthJointTaupVelSlip, "Eulerian_JointH_TaupVelSlip_vs_baralpp_%dX%dX%d"	,Filter,Filter,Filter);
				sprintf(charfilterwidthJointNVelSlip,	 "Eulerian_JointNreal_VelSlip_%dX%dX%d"			,Filter,Filter,Filter);
				sprintf(charfilterwidthJointNTaupVelSlip,"Eulerian_JointNreal_TaupVelSlip_%dX%dX%d"		,Filter,Filter,Filter);			

				outputfile = charfilterwidth;
				outputfileVelSlip = charfilterwidthVelSlip;	
				outputfileTaupVelSlip = charfilterwidthTaupVelSlip;				
				outputfileN = charfilterwidthN;		
				
				outputfileJointVelSlip = charfilterwidthJointVelSlip;			
				outputfileJointTaupVelSlip = charfilterwidthJointTaupVelSlip;				
				outputfileJointNVelSlip = charfilterwidthJointNVelSlip;			
				outputfileJointNTaupVelSlip = charfilterwidthJointNTaupVelSlip;	
			}

			else if (EulerianVslipBin_)
			{

				if ( nparticle_ < 0)
				{
					sprintf(charfilterwidth,                 "h_vs_baralpp_%dX%dX%d"                        	      ,Filter,Filter,Filter);
                        		sprintf(charfilterwidthVelSlip,          "velSlipEulerianManner_vs_baralpp_%dX%dX%d"                  ,Filter,Filter,Filter);
                        		sprintf(charfilterwidthTaupVelSlip,      "TaupVelSlipEulerianManner_vs_baralpp_%dX%dX%d"              ,Filter,Filter,Filter);
                       		 	sprintf(charfilterwidthN,                "Nrealization_%dX%dX%d"                        	      ,Filter,Filter,Filter);

                        		sprintf(charfilterwidthJointVelSlip,     "JointH_VelSlipEulerianManner_vs_baralpp_%dX%dX%d"           ,Filter,Filter,Filter);
                        		sprintf(charfilterwidthJointTaupVelSlip, "JointH_TaupVelSlipEulerianManner_vs_baralpp_%dX%dX%d"       ,Filter,Filter,Filter);
                        		sprintf(charfilterwidthJointNVelSlip,    "JointNreal_VelSlipEulerianManner_%dX%dX%d"                  ,Filter,Filter,Filter);
                        		sprintf(charfilterwidthJointNTaupVelSlip,"JointNreal_TaupVelSlipEulerianManner_%dX%dX%d"              ,Filter,Filter,Filter);
				}
				else
				{
					sprintf(charfilterwidth,		 "h_vs_baralpp_%dX%dX%d_npart=%d"					,Filter,Filter,Filter,nparticle_);
					sprintf(charfilterwidthVelSlip,		 "velSlipEulerianManner_vs_baralpp_%dX%dX%d_npart=%d"			,Filter,Filter,Filter,nparticle_);
					sprintf(charfilterwidthTaupVelSlip,	 "TaupVelSlipEulerianManner_vs_baralpp_%dX%dX%d_npart=%d"		,Filter,Filter,Filter,nparticle_);
					sprintf(charfilterwidthN,		 "Nrealization_%dX%dX%d_npart=%d"					,Filter,Filter,Filter,nparticle_);

					sprintf(charfilterwidthJointVelSlip, 	 "JointH_VelSlipEulerianManner_vs_baralpp_%dX%dX%d_npart=%d"		,Filter,Filter,Filter,nparticle_);
					sprintf(charfilterwidthJointTaupVelSlip, "JointH_TaupVelSlipEulerianManner_vs_baralpp_%dX%dX%d_npart=%d"	,Filter,Filter,Filter,nparticle_);
					sprintf(charfilterwidthJointNVelSlip,	 "JointNreal_VelSlipEulerianManner_%dX%dX%d_npart=%d"			,Filter,Filter,Filter,nparticle_);
					sprintf(charfilterwidthJointNTaupVelSlip,"JointNreal_TaupVelSlipEulerianManner_%dX%dX%d_npart=%d"		,Filter,Filter,Filter,nparticle_);
				}

				outputfile = charfilterwidth;
				outputfileVelSlip = charfilterwidthVelSlip;	
				outputfileTaupVelSlip = charfilterwidthTaupVelSlip;				
				outputfileN = charfilterwidthN;		
				
				outputfileJointVelSlip = charfilterwidthJointVelSlip;			
				outputfileJointTaupVelSlip = charfilterwidthJointTaupVelSlip;				
				outputfileJointNVelSlip = charfilterwidthJointNVelSlip;			
				outputfileJointNTaupVelSlip = charfilterwidthJointNTaupVelSlip;	
			}
			
			OFstream         str_gii(mesh.time().path()/outputRelativePath_/runTime.timeName()/outputfile);
			OFstream     str_velSlip(mesh.time().path()/outputRelativePath_/runTime.timeName()/outputfileVelSlip);
			OFstream str_TaupVelSlip(mesh.time().path()/outputRelativePath_/runTime.timeName()/outputfileTaupVelSlip);			
			OFstream       str_nreal(mesh.time().path()/outputRelativePath_/runTime.timeName()/outputfileN);

			OFstream                str_JointVelSlip(mesh.time().path()/outputRelativePath_/runTime.timeName()/outputfileJointVelSlip);
			OFstream            str_JointTaupVelSlip(mesh.time().path()/outputRelativePath_/runTime.timeName()/outputfileJointTaupVelSlip);			
			OFstream           str_JointNrealVelSlip(mesh.time().path()/outputRelativePath_/runTime.timeName()/outputfileJointNVelSlip);
			OFstream       str_JointNrealTaupVelSlip(mesh.time().path()/outputRelativePath_/runTime.timeName()/outputfileJointNTaupVelSlip);

			
			Pout << " 		   " << endl; 		
			Pout << " Writing the file " << outputfile << ", " << outputfileVelSlip      << ", " << outputfileTaupVelSlip      << ", " << outputfileN 
						     		   << ", " << outputfileJointVelSlip << ", " << outputfileJointTaupVelSlip << ", " << outputfileJointNVelSlip << ", " <<  outputfileJointNTaupVelSlip
			     << " into the folder: " << mesh.time().path()/outputRelativePath_/runTime.timeName() << endl;
			
			str_gii         		<< "# baralpp \t (1-H_x) \t (1-H_y) \t (1-H_z)	"     	<< nl;
			str_velSlip     		<< "# baralpp \t Vr_x \t Vr_y \t Vr_z		"  	<< nl;
			str_TaupVelSlip 		<< "# baralpp \t Tau_p*|Vr|  			"       << nl;
			str_nreal       		<< "# baralpp \t Nreal				"     	<< nl;
			
			str_JointVelSlip         	<< "# baralpp \t Vr_z 		\t (1-H_z)	"     	<< nl;
			str_JointTaupVelSlip     	<< "# baralpp \t Tau_p*|Vr| 	\t (1-H_z)	"  	<< nl;
			str_JointNrealVelSlip 		<< "# baralpp \t Vr_z		\t Nreal	"       << nl;
			str_JointNrealTaupVelSlip       << "# baralpp \t Tau_p*|Vr|	\t Nreal	"     	<< nl;

			for( int i = 0; i < nBin_; i++ )
			{
								
				str_nreal 	 << (i+1./2.)*(maxalpp_-minalpp_)/nBin_  << " " << numberCondalppFilter_[i][filterwidth] << nl;
				
				// Avoid to divide by zero
				if ( numberCondalppFilter_[i][filterwidth] == 0 ) numberCondalppFilter_[i][filterwidth] = 1; 
					  
				str_gii      	 << minalpp_+(i+1./2.)*(maxalpp_-minalpp_)/nBin_  
				          << " " << giiCondalppFilter_[i][0][filterwidth] / numberCondalppFilter_[i][filterwidth] 
					  << " " << giiCondalppFilter_[i][1][filterwidth] / numberCondalppFilter_[i][filterwidth] 
					  << " " << giiCondalppFilter_[i][2][filterwidth] / numberCondalppFilter_[i][filterwidth] 	<< nl;	
								       
				str_velSlip   	 << minalpp_+(i+1./2.)*(maxalpp_-minalpp_)/nBin_  
				          << " " << velSlipCondalppFilter_[i][0][filterwidth] / numberCondalppFilter_[i][filterwidth] 
					  << " " << velSlipCondalppFilter_[i][1][filterwidth] / numberCondalppFilter_[i][filterwidth] 
					  << " " << velSlipCondalppFilter_[i][2][filterwidth] / numberCondalppFilter_[i][filterwidth] 	<< nl;	

				str_TaupVelSlip  << minalpp_+(i+1./2.)*(maxalpp_-minalpp_)/nBin_  
				          << " " << taupVelSlipCondalppFilter_[i][filterwidth] / numberCondalppFilter_[i][filterwidth] 	<< nl;					
								       
								       
				if ( verbose_ )				       
				{
					Pout  	<< minalpp_+(i+1./2.)*(maxalpp_-minalpp_)/nBin_   
					 << " " << giiCondalppFilter_[i][0][filterwidth] 
					 << " " << giiCondalppFilter_[i][1][filterwidth]  
					 << " " << giiCondalppFilter_[i][2][filterwidth] 
					 << " " << numberCondalppFilter_[i][filterwidth] 						<< endl;
	
					Pout  	<< (minalpp_+i+1./2.)*(maxalpp_-minalpp_)/nBin_   
					 << " " << giiCondalppFilter_[i][0][filterwidth] / numberCondalppFilter_[i][filterwidth] 
					 << " " << giiCondalppFilter_[i][1][filterwidth] / numberCondalppFilter_[i][filterwidth] 
					 << " " << giiCondalppFilter_[i][2][filterwidth] / numberCondalppFilter_[i][filterwidth] 	<< endl;
				}
				
				for ( int j = 0; j < nBin_; j++)
				{
					
					str_JointNrealVelSlip 		<< minalpp_+(i+1./2.)*(maxalpp_-minalpp_)/nBin_  				<< " " 	
								  	<< minVelSlip_+(j+1./2.)*(maxVelSlip_-minVelSlip_)/nBin_			<< " "
					                        	<< numberVelSlipJointgiiCondalppFilter_[i][j][filterwidth]			<< nl; 
					
					str_JointNrealTaupVelSlip 	<< minalpp_+(i+1./2.)*(maxalpp_-minalpp_)/nBin_  				<< " " 	
								  	<< minTaupVelSlip_+(j+1./2.)*(maxTaupVelSlip_-minTaupVelSlip_)/nBin_  		<< " "
									<< numberTaupVelSlipJointgiiCondalppFilter_[i][j][filterwidth]  		<< nl;
									
					// Avoid to divide by zero
					if (     numberVelSlipJointgiiCondalppFilter_[i][j][filterwidth] == 0 )     numberVelSlipJointgiiCondalppFilter_[i][j][filterwidth] = 1; 				
					if ( numberTaupVelSlipJointgiiCondalppFilter_[i][j][filterwidth] == 0 ) numberTaupVelSlipJointgiiCondalppFilter_[i][j][filterwidth] = 1; 									

					str_JointVelSlip		<< minalpp_+(i+1./2.)*(maxalpp_-minalpp_)/nBin_  									<< " " 	
								  	<< minVelSlip_+(j+1./2.)*(maxVelSlip_-minVelSlip_)/nBin_  								<< " "
									<<    VelSlipJointgiiCondalppFilter_[i][j][filterwidth]/numberVelSlipJointgiiCondalppFilter_[i][j][filterwidth]  	<< nl;   

					str_JointTaupVelSlip		<< minalpp_+(i+1./2.)*(maxalpp_-minalpp_)/nBin_ 									<< " " 	
								  	<< minTaupVelSlip_+(j+1./2.)*(maxTaupVelSlip_-minTaupVelSlip_)/nBin_  							<< " "
									<< TaupVelSlipJointgiiCondalppFilter_[i][j][filterwidth]/numberTaupVelSlipJointgiiCondalppFilter_[i][j][filterwidth]  	<< nl; 				
											
				}
				
												

			} 
			
		}
		
		
	}
	
}	
