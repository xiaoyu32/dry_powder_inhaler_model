#include "coarsening.H"
#include "dragForce.H"
 
namespace Foam
{ 

 	void coarsening
	(
		cfdemCloud& parsm,
		const volScalarField& rho_,		
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volVectorField& tildeUf_,
		const volScalarField& p_,
		const volScalarField& barPg_,
		vectorField& pargii_,
		const labelListList& parList_,
		const bool& verbose_,
		const labelList& parIDInSubVolume_
	)
 	{
		// get viscosity field
		#ifdef comp
		    const volScalarField nufField = parsm.turbulence().mu()/rho_;
		#else
		    const volScalarField& nufField = parsm.turbulence().nu();
		#endif
		
		// Distance vector
		vector dist(0,0,0);
				
		// Parcel/Particle averaged force 
		vectorField parForce(parIDInSubVolume_.size(),vector(0,0,0));
				
		// Parcel/Particle averaged drag force
		vectorField parDragForce(parIDInSubVolume_.size(),vector(0,0,0));
		
		// Parcel/Particle averaged filtered relative velocity 
		vectorField parFilteredRelativeVelocity(parIDInSubVolume_.size(),vector(0,0,0));

		// Parcel/Particle averaged relative velocity
		vectorField parRelativeVelocity(parIDInSubVolume_.size(),vector(0,0,0));	
			
		volVectorField gradp_ = fvc::grad(p_);
		volVectorField gradbarPg_ = fvc::grad(barPg_);
			   
		interpolationCellPoint<vector> gradPgInterpolator_(gradp_);
		interpolationCellPoint<vector> gradbarPgInterpolator_(gradbarPg_);
		
		interpolationCellPoint<vector> UInterpolator_(Uf_);	
		interpolationCellPoint<vector> tildeUInterpolator_(tildeUf_);
		
		interpolationCellPoint<scalar> voidfractionInterpolator_(alpf_);
		
		vector gradPg_int(0,0,0);	
	   	vector gradbarPg_int(0,0,0);		
		vector position(0,0,0);
		vector Up(0,0,0);
		vector Ufluid(0,0,0);
		vector tildeUfluid(0,0,0);
		vector Ur(0,0,0);
		vector tildeUr(0,0,0);				
		scalar Volp(0);
		scalar ds(0);
		
		vector WenYuDrag(0,0,0);
		scalar voidfraction(0);
		scalar rhof(0);
		scalar nuf(0);
		
		label cellI(-1);
		
		
		// Number of parcell WARNING There is a problem here Why parcelIDInSubVolume_.size() != npLocalParcel
		int trimNpParcel = pow(10,10);		
		if( parIDInSubVolume_.size() < trimNpParcel ) trimNpParcel = parIDInSubVolume_.size();		
		if( parList_.size() < trimNpParcel ) trimNpParcel = parList_.size();	
																
		Pout << " Create parcel-> Trimming of max. number of parcel: " << trimNpParcel 
			 << " parcel sub-volume list: " << parIDInSubVolume_.size() << " parcel list size: " << parList_.size() << endl;					
		
		//for(int parI=0; parI < parIDInSubVolume_.size(); parI++ )
		for(int parI=0; parI < trimNpParcel; parI++ )
		//for(int parI=0; parI < parIDInSubVolume_.size(); parI++ )
		{ 						
			dist = vector(0,0,0);
			gradPg_int = vector(0,0,0);
			gradbarPg_int = vector(0,0,0);
			Ufluid = vector(0,0,0);
			tildeUfluid = vector(0,0,0);			
			Up = vector(0,0,0);
			Ur = vector(0,0,0);
			tildeUr = vector(0,0,0);			
			Volp = 0;
			WenYuDrag = vector(0,0,0);
			voidfraction = 0;
			nuf = 0;
			rhof = 0;	
						
			labelList parList_L = parList_[parI];	
			
			forAll(parList_L,particleI)
			{				
				label particleII = parList_L[particleI]; 
								
				cellI = parsm.cellIDs()[particleII][0];															
				
				if( cellI > -1 )
				{
					position = parsm.position(particleII); 
				
					gradPg_int = gradPgInterpolator_.interpolate(position,cellI); 
					gradbarPg_int = gradbarPgInterpolator_.interpolate(position,cellI); 
				
					Ufluid = UInterpolator_.interpolate(position,cellI); 
					tildeUfluid = tildeUInterpolator_.interpolate(position,cellI);	
				
					Up = parsm.velocity(particleII);
					Ur = Ufluid-Up;
					tildeUr = tildeUfluid-Up;
									 	
				    ds = 2*parsm.radius(particleII);
					Volp = ds*ds*ds*M_PI/6;
				
					// Calculate WenYu Drag 
					voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
					nuf = nufField[cellI];
					rhof = rho_[cellI];	
					WenYuDragForce(Ur,ds,rhof,nuf,voidfraction,WenYuDrag);	
							
					/*
					Pout << "parcel/particle " << parI << " by particles " << particleII << " cellI " << cellI 
						 << " voidfraction " << voidfraction << " Up " << Up << " WenYuDrag " << WenYuDrag << endl;
				
					*/				
																							
					// Parcel/Particle averaged resolved force
					for(int j=0;j<3;j++) 
					{		
													parForce[parI][j] += ( - Volp * gradPg_int[j] 
									                                       + Volp * gradbarPg_int[j]
																		   + WenYuDrag[j]	         )  ; 	
																		
								 parFilteredRelativeVelocity[parI][j] += tildeUr[j] ; 
				
								                parDragForce[parI][j] += WenYuDrag[j];
						
								         parRelativeVelocity[parI][j] += Ur[j] ; 						

					}			
				
					/*
					Pout << " force= " << parForce[parI] << " tildeUr= " << parFilteredRelativeVelocity[parI]
						 << " dragF= " << parDragForce[parI] << " Ur= " << parRelativeVelocity[parI] << endl;  		
					*/
				}
			}										
			
		}				

		// Parcel/Particle averaged filtered drag coefficient
		vectorField parbarBetai_(parIDInSubVolume_.size(),vector(0,0,0));
		
		// Drag correction with particle coarsening
		pargii_.resize(parIDInSubVolume_.size());
		
		//forAll(parList_,parII)
		for(int ii=0; ii < parIDInSubVolume_.size(); ii++ )
		{ 	
			for(int j=0;j<3;j++) 
			{				
			 parbarBetai_[ii][j] =        parForce[ii][j] / ( parFilteredRelativeVelocity[ii][j] + SMALL ) ;
			      pargii_[ii][j] =    parbarBetai_[ii][j] *           parRelativeVelocity[ii][j] 
				                       / ( parDragForce[ii][j] + SMALL ) ;				
			}
		}		
		
		if(verbose_)
		{
			int index = 0;			
			Pout << " Parcel/Particle barBetai " << parbarBetai_[index] << endl;
			Pout << " Parcel/Particle gii      " << pargii_[index]      << endl;		
		}
	
	}
	
}	
