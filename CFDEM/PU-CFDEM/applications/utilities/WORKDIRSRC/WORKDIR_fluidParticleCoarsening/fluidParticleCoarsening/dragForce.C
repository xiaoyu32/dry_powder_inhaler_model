#include "dragForce.H"

namespace Foam
{	
	void WenYuDragForce
	(		
		const vector& Ur_,
		const scalar& ds_,
		const scalar& rho_,
		const scalar& nuf_,
		const scalar& voidfraction_,
		vector& WenYuDrag_
	)
	{
		scalar Rep(0);
		scalar magUr(0);
		scalar CD(0);
		scalar betaP(0);
		
		scalar Volp(0);	
		
		magUr = mag(Ur_);		
		Volp = ds_*ds_*ds_*M_PI/6;

        if(ds_<=0) Pout << " ds " << ds_ << endl;
		
		if (magUr > 0 && ds_ > 0)
		{
			// calc particle Re number
			Rep = voidfraction_*ds_*magUr/(nuf_+SMALL);
            
			// calc CD
			if (Rep < 1000)
			{
				CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction_,-2.65);
			}
			else
			{
				CD = 0.44*pow(voidfraction_,-2.65); 
			}  	

			// calc drag coefficient 
			betaP = 3./4.*rho_*CD*magUr/ds_;  	

			// calc particle's drag
			WenYuDrag_ = Volp*betaP*Ur_;

		}else
		{
			WenYuDrag_ = vector(0,0,0);
		}	
						
	}

	void CalculateDragForce
	(
		cfdemCloud& sm,
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volScalarField& rho_,
		const bool& verbose_,
		vectorField& DragForce_,
		const labelList& partIDInSubVolume_
	)
	{
		
		// get viscosity field
		#ifdef comp
		    const volScalarField nufField = sm.turbulence().mu()/rho_;
		#else
		    const volScalarField& nufField = sm.turbulence().nu();
		#endif

		// Local variables	
		label  cellI(-1);
		vector drag(0,0,0);
		vector Ufluid(0,0,0);
		
		vector position(0,0,0);
		scalar voidfraction(1);
		
		vector Up(0,0,0);
		vector Ur(0,0,0);
		scalar ds(0);
		
		scalar nuf(0);
		scalar rhof(0);
		
		vector WenYuDrag(0,0,0);
		
		interpolationCellPoint<scalar> voidfractionInterpolator_(alpf_);
		interpolationCellPoint<vector> UInterpolator_(Uf_);	
				
		// 
		//_AO_Parallel
		//DragForce_.resize(sm.numberOfParticles());
		DragForce_.resize(partIDInSubVolume_.size());
					
		//for(int index = 0; index <  sm.numberOfParticles(); index++)
	    for(int ii =0; ii < partIDInSubVolume_.size(); ii++)
		{
			int index = partIDInSubVolume_[ii];
			cellI = sm.cellIDs()[index][0];
			drag = vector(0,0,0);
			Ufluid = vector(0,0,0);
			WenYuDrag = vector(0,0,0);
			//DragForce_[index] = vector(0,0,0);
			DragForce_[ii] = vector(0,0,0);
			    
			if (cellI > -1) // particle Found
			{
				position = sm.position(index);
			
				if ( alpf_[cellI] > 1. ) Pout << " voidfraction > 1 " << alpf_[cellI] << endl;
			
				voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
				Ufluid = UInterpolator_.interpolate(position,cellI);
			
				if ( voidfraction > 1. ) 
				{
						Pout << " Int. voidfraction > 1 " << " value= " << voidfraction;
						voidfraction = alpf_[cellI];
						Pout << " mod. value = " << voidfraction << endl;
				}		
				Up = sm.velocity(index);
			
				Ur = Ufluid-Up;				
				ds = 2*sm.radius(index);		
				rhof = rho_[cellI];
				nuf = nufField[cellI];
			
				// Drag force
				WenYuDragForce(Ur,ds,rhof,nuf,voidfraction,WenYuDrag);
						
				if(verbose_ && index <= 1)
				{
					Info << "" << endl;
					Pout << " index = " << index << endl;
					Pout << " position = " << position << endl; 
					Pout << " Up = " << Up << endl;
					Pout << " Ur = " << Ur << endl;
					Pout << " dp = " << ds << endl;
					Pout << " rho = " << rhof << endl;
					Pout << " nuf = " << nuf << endl;
					Pout << " voidfraction = " << voidfraction << endl;
					Pout << " drag = " << WenYuDrag << endl;
					Info << " " << endl;
				}
			}	
			
			//for(int j=0;j<3;j++) DragForce_[index][j] = drag[j];		
			for(int j=0;j<3;j++) DragForce_[ii][j] = WenYuDrag[j];	
		}	
	}
	
	void EulerianParticleVelocityForce
	(
		cfdemCloud& sm,			
		const fvMesh& mesh,
		volVectorField&	Up_,
		const vectorField& DragForce_,
		volVectorField& MappedDragForce_,
		const labelList& partIDInSubVolume_		
	)
	{		
		// Neighbouring cells
		CPCCellToCellStencil neighbourCells(mesh);	
		
		label cellID;
		vector position(0,0,0);
		vector particleVel(0,0,0);
		vector particleDrag(0,0,0);
		scalar dist_s(0);
		scalar sumWeights(0);
		
		scalarField               weightScalar(27,scalar(0.0));
		//Field <Field <scalar> >   particleWeights(sm.numberOfParticles(),weightScalar);		
		Field <Field <scalar> >   particleWeights(partIDInSubVolume_.size(),weightScalar);
		
		
		//for(int index = 0; index <  sm.numberOfParticles(); index++)
	    	for(int ii =0; ii < partIDInSubVolume_.size(); ii++)
		{
			int index = partIDInSubVolume_[ii];		
			cellID = sm.cellIDs()[index][0];
			position = sm.position(index);			    
			particleVel  = sm.velocity(index);
			particleDrag = DragForce_[index];
			
			//if (cellID > -1)  // particle centre is in domain
            //		{
				labelList& cellsNeigh = neighbourCells[cellID];
				sumWeights = 0;
				dist_s = 0;
								
				forAll(cellsNeigh,ii)
				{
					// Find distances between particle and neighbouring cells					
					dist_s = mag(sm.mesh().C()[cellsNeigh[ii]]-position)/pow(sm.mesh().V()[cellsNeigh[ii]],1./3.);
															
					if(dist_s <= 0.5)
					{		
						particleWeights[index][ii] =  1./4.*pow(dist_s,4)-5./8.*pow(dist_s,2)+115./192.;
					}
					else if (dist_s > 0.5 && dist_s <= 1.5)
					{		
						particleWeights[index][ii] = -1./6.*pow(dist_s,4)+5./6.*pow(dist_s,3)-5./4.*pow(dist_s,2)+5./24.*dist_s+55./96.;
					}
					else if (dist_s > 1.5 && dist_s <= 2.5)
					{		
						particleWeights[index][ii] =  pow(2.5-dist_s,4)/24.;
					}
					else
					{		
						particleWeights[index][ii] = 0;
					}
					
					sumWeights += particleWeights[index][ii];
						
				}	
				
				forAll(cellsNeigh,ii)
				{	
					if ( sumWeights != 0 )
					{
						Up_[cellID] 	         +=  particleVel*particleWeights[index][ii]/sumWeights;
						MappedDragForce_[cellID] += particleDrag*particleWeights[index][ii]/sumWeights;
					}
					else
					{
						Up_[cellID] 		 = vector(0,0,0);
						MappedDragForce_[cellID] = vector(0,0,0);	
					}
				}	
			//}	
		}
			
									
	}
	void CalculateEulerianDragForce
	(		
		cfdemCloud& sm,
		const fvMesh& mesh,
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volVectorField& Up_,
		const volScalarField& rho_,
		volVectorField& EulerianDragForce_	
	)
	{
		// get viscosity field
		#ifdef comp
		    const volScalarField nufField = sm.turbulence().mu()/rho_;
		#else
		    const volScalarField& nufField = sm.turbulence().nu();
		#endif
		
		scalar ds(0);
		scalar rhof(0);
		scalar nuf(0);
		scalar voidfraction(1);
		scalar alps(0);
		vector Ur(0,0,0);
		vector WenYuDrag(0,0,0);
		
		int partID = 0; // Monodisperse case, be careful if for the cases with different particle diameters
		forAll(mesh.cells(),cellI)
		{
			Ur = Uf_[cellI]-Up_[cellI];
			ds = 2*sm.radius(partID);
			rhof = rho_[cellI];		
			nuf = nufField[cellI];
			voidfraction = alpf_[cellI];
			
			// Drag force
			WenYuDragForce(Ur,ds,rhof,nuf,voidfraction,WenYuDrag);
			
			alps = 1-voidfraction+SMALL;
			
			EulerianDragForce_[cellI] = alps*WenYuDrag;

		}		
	}




	void filteredEulerEulerDragCoefficient
	(
		cfdemCloud& sm,
		const bool& verbose_,
		const fvMesh& mesh,
		const labelListListList& stencillistfilter,
		const int& filterwidth,		
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volVectorField& Up_,		
		const volScalarField& baralpf_,
		const volVectorField& tildeUf_,
		const volVectorField& tildeUp_,		
		const volScalarField& barPg_,
		const volScalarField& p_,
		const volScalarField& rho_,
		const volVectorField& DragForce_,
		volVectorField& ResolvedDragForce_,
		vectorField& gii_
	)

	{
		
		volVectorField Beta_fil
		(   
		    IOobject
		    (
        		"Beta_fil",
        		mesh.time().timeName(),
        		mesh,
        		IOobject::NO_READ,
        		IOobject::NO_WRITE
		    ),
		    mesh,
		    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
		); 	
				
		Beta_fil  = ( -    rho_*(1.-alpf_)*fvc::grad(p_)
			      + rho_*(1.-baralpf_)*fvc::grad(barPg_)
			      + DragForce_   		       	     ); 	      	
		
		// Correction in vector array
		forAll(mesh.cells(),ii)
		{

			
			if( ( 1.-alpf_[ii] ) > SMALL ) // If there is no particle in the cell then no drag & correction ...   
			{
				Beta_fil[ii].component(vector::X) /= ( tildeUf_[ii].component(vector::X) - tildeUp_[ii].component(vector::X) );
				Beta_fil[ii].component(vector::Y) /= ( tildeUf_[ii].component(vector::Y) - tildeUp_[ii].component(vector::Y) );
				Beta_fil[ii].component(vector::Z) /= ( tildeUf_[ii].component(vector::Z) - tildeUp_[ii].component(vector::Z) );
				
				gii_[ii][0] = Beta_fil[ii].component(vector::X) * ( Uf_[ii].component(vector::X) - Up_[ii].component(vector::X) ) / DragForce_[ii].component(vector::X);
				gii_[ii][1] = Beta_fil[ii].component(vector::Y) * ( Uf_[ii].component(vector::Y) - Up_[ii].component(vector::Y) ) / DragForce_[ii].component(vector::Y);
				gii_[ii][2] = Beta_fil[ii].component(vector::Z) * ( Uf_[ii].component(vector::Z) - Up_[ii].component(vector::Z) ) / DragForce_[ii].component(vector::Z);
			}
			else	
			{
				gii_[ii][0] = 0;
				gii_[ii][1] = 0;
				gii_[ii][2] = 0;
			}
		}
		
		/*
		forAll(mesh.cells(),ii)
		{
				Info << " Diff resolved press contr. = " << - rho_[ii]*alpf_[]*fvc::grad(p_)[ii].component(vector::X) + rho_*baralpf_*fvc::grad(barPg_) 
		}
		*/
		
													 
   	}
}	
