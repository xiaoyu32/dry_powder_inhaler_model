#include "applyBin.H"

namespace Foam
{
	void applyBins
	(
		cfdemCloud& sm,
		const fvMesh& mesh,
		const int& filterwidth,
		const int& maxfilterwidth,
		const int& nBin_,
		const scalar& minalpp_,
		const scalar& maxalpp_,
		const scalar& minVelSlip_,
		const scalar& maxVelSlip_,
		const scalar& minTaupVelSlip_,
		const scalar& maxTaupVelSlip_,		
		const scalar& rhop,
		const volScalarField& rho_,
		const volScalarField& alpf_,		
		const volScalarField& baralpf_,
		const volVectorField& Uf_,		
		const volVectorField& tildeUf_,
		const volVectorField& Up_,		
		const volVectorField& tildeUp_,		
		vectorField& gii_,
		Field< Field < Field <scalar> > >& 	giiCondalppFilter_,
		Field< Field < Field <scalar> > >& 	velSlipCondalppFilter_,
		Field < Field <scalar> >& 		taupVelSlipCondalppFilter_,
		Field < Field <scalar> >& 		numberCondalppFilter_, 		
		Field< Field < Field <scalar> > >&	VelSlipJointgiiCondalppFilter_,
		Field< Field < Field <scalar> > >&	TaupVelSlipJointgiiCondalppFilter_,	
		Field <Field < Field <scalar> > >& 	numberVelSlipJointgiiCondalppFilter_,
		Field <Field < Field <scalar> > >& 	numberTaupVelSlipJointgiiCondalppFilter_,
		bool EulerEulerFiltering_,
		bool EuerianVslipbin_,
		const labelList& partIDInSubVolume_ 
	)
	{			
		// get viscosity field
		#ifdef comp
		    const volScalarField nufField = sm.turbulence().mu()/rho_;
		#else
		    const volScalarField& nufField = sm.turbulence().nu();
		#endif
		
		//const scalar minalpp_ = 0.;
		if ( maxalpp_ == 0)
		{
			FatalError << " Max. solid volume fraction equal to = " << maxalpp_ << abort(FatalError);
		}	
		
		// Binning particle solid volume fraction
		const scalar binDx_   = (maxalpp_-minalpp_)/nBin_;	
		// Binning relative velocity 
		const scalar binVelSlip_ = (maxVelSlip_-minVelSlip_)/nBin_;
		// Binning tau_p |Vr|
		const scalar binTaupVelSlip_= (maxTaupVelSlip_-minTaupVelSlip_)/nBin_;
		
		label cellI(-1);
		vector position(0,0,0);
		vector intertildeUf_(0,0,0);
		vector Up(0,0,0);
		vector Ur(0,0,0);
		scalar baralpp(0);
		scalar CD(0);		
		scalar nuf(0);
		scalar rhof(0);
		scalar baralpf(0);
		scalar magUr(0);
		scalar ds(0);
		scalar Rep(0);
		scalar taup(0);		
		
		int maxIndex;		
		
		interpolationCellPoint<scalar> baralpfInterpolator_(baralpf_);		
		interpolationCellPoint<vector> tildeUInterpolator_(tildeUf_);
		interpolationCellPoint<vector> tildeUpInterpolator_(tildeUp_);
				
		if(EulerEulerFiltering_) // Lagrangian filtering
		{
			maxIndex = mesh.cells().size();
		}
		//else	maxIndex = sm.numberOfParticles();
		else maxIndex = partIDInSubVolume_.size(); 	
									
		//for(int index = 0; index <  sm.numberOfParticles(); index++)
		for(int ii =0; ii < maxIndex; ii++)
		{
			int index = partIDInSubVolume_[ii];			
			
			if(EulerEulerFiltering_) 
			{
				cellI = ii; //index;
			}
			else cellI = sm.cellIDs()[index][0];
			
			//if (cellI > -1) // particle Found
			//{

				if(EulerEulerFiltering_)
				{
					Up = Up_[cellI];
					position = mesh.C()[cellI];	
				}
				else if(EuerianVslipbin_)
				{
					Up = Up_[cellI];
					position      = sm.position(index);
				}
				else
				{
					Up            = sm.velocity(index);
					position      = sm.position(index);				
				}				
				
				baralpf       =    baralpfInterpolator_.interpolate(position,cellI) ;
				baralpp       = 1.-baralpfInterpolator_.interpolate(position,cellI) ;
				intertildeUf_ =     tildeUInterpolator_.interpolate(position,cellI) ;
				
				Ur = intertildeUf_-Up;
				
				if(EuerianVslipbin_) Ur = intertildeUf_ - tildeUpInterpolator_.interpolate(position,cellI);
				
				/*
				if(EulerEulerFiltering_)
				{
					ds = 2*sm.radius(index);		
				} else	ds = 2*sm.radius(0); // Be careful for monodisperse
				*/
				ds = 2*sm.radius(0); // Be careful for monodisperse
				
				nuf = nufField[cellI];
				rhof = rho_[cellI];
				magUr = mag(Ur);

				if (magUr > 0)
				{
					// calc particle Re number
					Rep = baralpf*ds*magUr/(nuf+SMALL);

					// calc CD
					if (Rep < 1000)
					{
						CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(baralpf,-2.65);
					}
					else
					{
						CD = 0.44*pow(baralpf,-2.65); 
					}  
				}

				// Relaxation time
				taup = 4./3.*rhop/rhof*ds/CD;
				
				label binI   = floor((baralpp-minalpp_)/binDx_);
								
                if ( binI < nBin_ && binI > -1 ) // Avoid to sample out of the defined binning interval
				{					
					for(int j=0;j<3;j++) 
					{	
					      //giiCondalppFilter_[binI][j][filterwidth] += gii_[index][j];
					      giiCondalppFilter_[binI][j][filterwidth] += gii_[ii][j];				          
					  velSlipCondalppFilter_[binI][j][filterwidth] += intertildeUf_[j]-Up[j];  										
					}
				         taupVelSlipCondalppFilter_[binI][filterwidth] += taup;				      					      					      
					      numberCondalppFilter_[binI][filterwidth] += 1;
				}
				
				// Along z-direction
				label binVelSlipI     = floor((intertildeUf_[2]-Up[2]-minVelSlip_)/binVelSlip_); 

				label binTaupVelSlipI = floor((taup*magUr-minTaupVelSlip_)/binTaupVelSlip_);

				if ( binI < nBin_ && binVelSlipI < nBin_ && binVelSlipI > -1 && binI > -1 )
				{
				   	         //VelSlipJointgiiCondalppFilter_[binI][binVelSlipI][filterwidth] += gii_[index][2];
				   	         VelSlipJointgiiCondalppFilter_[binI][binVelSlipI][filterwidth] += gii_[ii][2];
					   numberVelSlipJointgiiCondalppFilter_[binI][binVelSlipI][filterwidth] += 1;
				} 

				if ( binI < nBin_ && binTaupVelSlipI < nBin_ && binTaupVelSlipI > -1 && binI > -1 )
				{    
					 //TaupVelSlipJointgiiCondalppFilter_[binI][binTaupVelSlipI][filterwidth] += gii_[index][2];
					 TaupVelSlipJointgiiCondalppFilter_[binI][binTaupVelSlipI][filterwidth] += gii_[ii][2];
				   numberTaupVelSlipJointgiiCondalppFilter_[binI][binTaupVelSlipI][filterwidth] += 1;
				}


			//}
		}

			
	}
}	
