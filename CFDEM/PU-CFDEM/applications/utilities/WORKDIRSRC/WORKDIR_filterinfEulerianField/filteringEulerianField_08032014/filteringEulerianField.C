/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    postCalc

Description
    Generic wrapper for calculating a quantity at each time

\*---------------------------------------------------------------------------*/

#include "filteringEulerianField.H"
#include "timeSelector.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

	void createStencils
	(
		const fvMesh& mesh,
		const int& filterwidth_,
		const labelListList& mA_,
		const int& max_Ax_,
		const int& max_Ay_,
		const int& max_Az_,
		labelListList& stencillist_
	)
	{
		labelList nStencils(pow(2*(filterwidth_+1)+1,3),-1);
		stencillist_= nStencils;
		
		label cellII = 0;
		
		forAll(mesh.cells(),cellI)
		{		
			int count = 0;	
			//Info << " Cell = " << cellI << endl;
			for ( label kwidth = -(filterwidth_+1); kwidth <= filterwidth_+1; kwidth++)
			{
				for ( label jwidth = -(filterwidth_+1); jwidth <= filterwidth_+1; jwidth++)
				{
					for ( label iwidth = -(filterwidth_+1); iwidth <= filterwidth_+1; iwidth++)
					{
						/*
						Info << mA_[cellI]      << " " << (mA_[cellI][0]+iwidth+max_Ax_+1) % (max_Ax_+1)
									<< " " << (mA_[cellI][1]+jwidth+max_Ay_+1) % (max_Ay_+1)
									<< " " << (mA_[cellI][2]+kwidth+max_Az_+1) % (max_Az_+1)
									<< " " << endl;
						
						*/
						cellII =  ((mA_[cellI][0]+iwidth+max_Ax_+1) % (max_Ax_+1))
						         +((mA_[cellI][1]+jwidth+max_Ay_+1) % (max_Ay_+1))*(max_Ax_+1)
							 +((mA_[cellI][2]+kwidth+max_Az_+1) % (max_Az_+1))*(max_Ax_+1)*(max_Ay_+1) ;
						
						/*
						Info << " Cell II " << cellII << endl;
						*/
													  			
						stencillist_[cellI][count] = cellII;					
						count++;
						
					}	
				}
			}
		}
	
	}

    	void constructfilter
        (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh,
                const int& minfilterwidth,
                const int& maxfilterwidth,
                const int& filterincrement,
                labelListListList& stencillistfilter,
                const fileName& outputRelativePath_
        )
	
	{
	
		const pointField& pp = mesh.points();

		// Min, max x-coordinates	
		scalar min_x = Foam::min(pp & vector(1,0,0));
		scalar max_x = Foam::max(pp & vector(1,0,0));

		// Min, max y-coordinates		
		scalar min_y = Foam::min(pp & vector(0,1,0));
		scalar max_y = Foam::max(pp & vector(0,1,0));

		// Min, max z-coordinates		
		scalar min_z = Foam::min(pp & vector(0,0,1));
		scalar max_z = Foam::max(pp & vector(0,0,1));	

		Info << "Domain " << " x=[" << min_x << ":" << max_x << "]" << " y=[" << min_y << ":" << max_y << "]" << " z=[" << min_z << ":" << max_z << "]" << endl;

		const scalar cell_volume_ref = mesh.V()[0];
		forAll(mesh.cells(),cellI)
		{
			if ( abs(cell_volume_ref-mesh.V()[cellI]) > SMALL)  FatalError<< "Non-uniform mesh !!! " << abort(FatalError);
		}

		const scalar delta = pow(cell_volume_ref,1./3);
		Info << "Uniform mesh,  volume of a cell = " << cell_volume_ref << endl;
		Info << "Delta_x = Delta_y = Delta_z = " << delta << endl;


		// Boundary patches
		// Find coupled patches
		const polyBoundaryMesh& patches = mesh.boundaryMesh();

		label nCoupledPatches = 0; 	
		forAll(patches, patchI)
		{
			const polyPatch& pp = patches[patchI];

			if (pp.coupled())
			{
			    nCoupledPatches++;
			}					    
		}
		
		if ( nCoupledPatches < 6 ) FatalError<< " The domain is not fully periodic !!! " << abort(FatalError);
		Info << "Number of coupled patches = " << nCoupledPatches << endl;
		Info << " " << endl;
		
		// Create cell matrix
		labelList mACell(3,-1);
		labelListList  mA(mesh.nCells(),mACell);
		int mAx(0); int mAy(0); int mAz(0);
		int max_mAx(0); int max_mAy(0); int max_mAz(0);		
		forAll(mesh.cells(),cellI)
		{
			mAx = floor((mesh.C()[cellI][0]-min_x)/delta);
			mAy = floor((mesh.C()[cellI][1]-min_y)/delta);
			mAz = floor((mesh.C()[cellI][2]-min_z)/delta);			
			mA[cellI][0] = mAx;
			mA[cellI][1] = mAy;
			mA[cellI][2] = mAz;

			if ( mAx > max_mAx ) max_mAx = mAx; 
			if ( mAy > max_mAy ) max_mAy = mAy; 
			if ( mAz > max_mAz ) max_mAz = mAz;			
		}	

		//Info << max_mAx << " " << max_mAy << " " << max_mAz << endl;
		
		char charfilterwidth[100];
		
		// Find stencils						
		labelListList stencillist(mesh.cells().size());
		if ( maxfilterwidth !=0 ) 
		{
			stencillistfilter.resize(maxfilterwidth*stencillist.size());				
		}
		else // avoid list with zero length for filter 3X3X3
		{
			stencillistfilter.resize(stencillist.size());	
		}
		// Create stencil list folder if not exist
		if( !isDir(mesh.time().path()/outputRelativePath_) )
		{
			mkDir(mesh.time().path()/outputRelativePath_);													
 
			for( int filterwidth = minfilterwidth; filterwidth <= maxfilterwidth; filterwidth += filterincrement )
			{
				int Filter = 2*(filterwidth+1)+1;
				sprintf(charfilterwidth, "stencils_%dX%dX%d",Filter,Filter,Filter);
				fileName outputfile(charfilterwidth);		

				if ( !isFile(mesh.time().path()/outputRelativePath_/outputfile) )
				{
					Info << "Creating stencil list for " << Filter << "X" << Filter << "X" << Filter << " filter " << endl;
					OFstream str_stencil(mesh.time().path()/outputRelativePath_/outputfile);
				
					// Call multipleCelltoPoints
					createStencils
					(
                				mesh,
						filterwidth,
						mA,
						max_mAx,
						max_mAy,
						max_mAz,
						stencillist								
					);

					forAll(mesh.cells(),cellI)
					{										
						labelList cellSt = stencillist[cellI];
						str_stencil << cellSt.size() << "( " ;

						forAll(cellSt,StI)
						{
							str_stencil << cellSt[StI] << " " ;
						}
						str_stencil << ")" << nl;
					}

					stencillistfilter[filterwidth] = stencillist;
				}	

			}	
		}
		else
		{
			for( int filterwidth = minfilterwidth; filterwidth <= maxfilterwidth; filterwidth += filterincrement )
			{
				int Filter = 2*(filterwidth+1)+1;
				sprintf(charfilterwidth, "stencils_%dX%dX%d",Filter,Filter,Filter);
				fileName outputfile(charfilterwidth);
				
				
				if ( isFile(mesh.time().path()/outputRelativePath_/outputfile) )
				{					
					Info << "Reading stencil list for " << Filter << "X" << Filter << "X" << Filter << " filter " << endl;	
					IFstream str_stencil(mesh.time().path()/outputRelativePath_/outputfile);
					forAll(mesh.cells(),cellI)
					{
						str_stencil >> stencillist[cellI];
					}
					
					stencillistfilter[filterwidth] = stencillist;																			
				}
				else
				{
					Info << "Creating stencil list for " << Filter << "X" << Filter << "X" << Filter << " filter " << endl;
					OFstream str_stencil(mesh.time().path()/outputRelativePath_/outputfile);
					
					// Call multipleCelltoPoints
					createStencils
					(
                				mesh,
						filterwidth,
						mA,
						max_mAx,
						max_mAy,
						max_mAz,
						stencillist		
					);
					
					forAll(mesh.cells(),cellI)
					{
										
						labelList cellSt = stencillist[cellI];
						str_stencil << cellSt.size() << "( " ;
						
						forAll(cellSt,StI)
						{
							str_stencil << cellSt[StI] << " " ;
						}
						str_stencil << ")" << nl;

					}
					
					stencillistfilter[filterwidth] = stencillist;					
				}				
												
			}			
		}		

	}
	
	void readEulerianVariables
	(
	    	const argList& args, 
	    	const Time& runTime, 
	    	const fvMesh& mesh,
		volScalarField& alpf_,
		volVectorField& Uf_,
		volScalarField& rho_,
		volScalarField& p_	
	)
	{
		// Read gas density
		IOobject rhoheader
		(
		   "rho",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		   
			Info<< " 	Reading rho" << endl;
			volScalarField density_(rhoheader,mesh);	
			
			rho_ = density_ ;		

		// Read volume fraction of gas
		IOobject voidfractionheader
		(
		   "voidfraction",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);

			Info<< " 	Reading voidfraction" << endl;
			volScalarField voidfraction_(voidfractionheader,mesh);	
			alpf_ = voidfraction_ ;		

		// Read gas velocity 
		IOobject Uheader
		(
		   "U",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		   
			Info<< " 	Reading U" << endl;
			volVectorField U_(Uheader,mesh);
			Uf_ = U_ ;
			
		
		// Read gas pressure
		IOobject pheader
		(
		   "p",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
	   
			Info<< " 	Reading Pg" << endl;
			volScalarField Pg_(pheader,mesh);
			p_ = Pg_ ;	
	}

	void EulerianParticleVelocity
	(
		cfdemCloud& sm,			
		const fvMesh& mesh,
		volVectorField&	Up_			
	)
	{
		scalarField npartCell(mesh.cells().size(),scalar(0));	 
		label cellI = -1;
		for(int index = 0; index <  sm.numberOfParticles(); index++)
	    	{
			cellI = sm.cellIDs()[index][0];
			    
			if (cellI > -1) // particle Found
			{
				Up_[cellI] += sm.velocity(index);
				npartCell[cellI] ++;
			}
		}
		
		forAll(mesh.cells(),cellI)
		{
			if ( npartCell[cellI] != 0 )
			{
				Up_[cellI] /= npartCell[cellI];
			}
			else
			{
				Up_[cellI] = vector(0,0,0);	
			}
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
		scalar rho(0);
		scalar nuf(0);
		scalar Rep(0);
		scalar alps(0);
		scalar voidfraction(1);
		scalar magUr(0);
		scalar CD(0);
		vector Ur(0,0,0);
		
		int partID = 0; // Monodisperse case, be careful if for the cases with different particle diameters
		forAll(mesh.cells(),cellI)
		{
			Ur = Uf_[cellI]-Up_[cellI];
			ds = 2*sm.radius(partID);		
			nuf = nufField[cellI];
			rho = rho_[cellI];
			magUr = mag(Ur);
			voidfraction = alpf_[cellI];
			alps = 1-voidfraction+SMALL;

			if (magUr > 0)
			{
				// calc particle Re number
				Rep = voidfraction*ds*magUr/(nuf+SMALL);

				// calc CD
				if (Rep < 1000)
				{
					CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-2.65);
				}
				else
				{
					CD = 0.44*pow(voidfraction,-2.65); 
				}  	

				// calc particle's drag
				EulerianDragForce_[cellI] = 3./4.*alps*rho*CD*magUr/ds*Ur;

			}

		}		
	}
	    
	void filteringEulerianVariables
	(
		const argList& args, 
		const Time& runTime, 
		const fvMesh& mesh,
		const labelListListList& stencillistfilter,
		const int& filterwidth,
		const volScalarField& voidfraction_,
		const volVectorField& U_,
		const volScalarField& p_,
		volScalarField& baralpf_,
		volVectorField& tildeUf_,
		volScalarField& barPg_,
		const volVectorField& Up_,
		volVectorField& tildeUp_,
		const bool EulerianParticleVelocity_  
	)
	{

		char charfilterwidth[100]; 
		int Filter = 2*(filterwidth+1)+1;
		
		// Filtering volume fraction of gas phase		
		Info<< " 	Filtering voidfraction" << endl;
		sprintf(charfilterwidth, "barvoidfraction_%dX%dX%d",Filter,Filter,Filter);	

		volScalarField barvoidfraction
		(
		    IOobject
		    (
		        charfilterwidth,
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ
		    ),
			voidfraction_	
		);
				
		
		forAll(mesh.cells(),cellI)
		{
		    scalar total_volume = 0;
			barvoidfraction[cellI] = 0; 
			const labelList fcellI = stencillistfilter[filterwidth][cellI]; 
											
			forAll(fcellI,filtercellI)
			{
				total_volume           +=       mesh.V()[fcellI[filtercellI]];
				barvoidfraction[cellI] +=  voidfraction_[fcellI[filtercellI]] 
					                          * mesh.V()[fcellI[filtercellI]]; 
			}
		        if( total_volume > 0 )
			{
				barvoidfraction[cellI] /= total_volume; 
			}
			else 
			{
				barvoidfraction[cellI] = 0;
			}
		} 
		
		Info<< " 	Writing filtered voidfraction" << endl;
	    	barvoidfraction.write();		
		baralpf_ = barvoidfraction;
		
		// Filtering gas velocity		
		Info<< " 	Filtering U" << endl;
		sprintf(charfilterwidth, "tildeU_%dX%dX%d",Filter,Filter,Filter);	
		
		volVectorField tildeU
		(
		    IOobject
		    (
		        charfilterwidth,
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ
		    ),
		    U_
		);
	
		forAll(mesh.cells(),cellI)
		{
		    	scalar total_volume = 0;
			tildeU[cellI] = vector(0,0,0); 			
			const labelList fcellI = stencillistfilter[filterwidth][cellI]; 
											
			forAll(fcellI,filtercellI)
			{
				total_volume  +=       mesh.V()[fcellI[filtercellI]];
				tildeU[cellI] +=  voidfraction_[fcellI[filtercellI]] 
					            *            U_[fcellI[filtercellI]] 
								*      mesh.V()[fcellI[filtercellI]]; 
			}
	        if( total_volume > 0 )
			{
				tildeU[cellI] /= total_volume; 
			
				if ( barvoidfraction[cellI] > 0 )
				{
					tildeU[cellI] /= barvoidfraction[cellI];
				}
			}
			else 
			{
				tildeU[cellI] = vector(0,0,0);
			}
		} 
				
		Info<< " 	Writing filtered U" << endl;
	    	tildeU.write();		
		tildeUf_ = tildeU;
		
		// Filtering volume fraction of gas phase		
		Info<< " 	Filtering gas pressure" << endl;
		sprintf(charfilterwidth, "barPg_%dX%dX%d",Filter,Filter,Filter);	
		
		volScalarField barp
		(
		    IOobject
		    (
		        charfilterwidth,
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ
		    ),
		    p_	
		);
		
		forAll(mesh.cells(),cellI)
		{
		    	scalar total_volume = 0;
			barp[cellI] = 0; 
			const labelList fcellI = stencillistfilter[filterwidth][cellI]; 
											
			forAll(fcellI,filtercellI)
			{
				total_volume +=   mesh.V()[fcellI[filtercellI]];
				barp[cellI]  +=         p_[fcellI[filtercellI]] 
					            * mesh.V()[fcellI[filtercellI]]; 
			}
		        if( total_volume > 0 )
			{
				barp[cellI] /= total_volume; 
			}
			else 
			{
				barp[cellI] = 0;
			}
		} 
		
		Info<< " 	Writing filtered gas pressure" << endl;
	    	barp.write();		
		barPg_ = barp;	
		
		
		
		if(EulerianParticleVelocity_)
		{			
			// Filtering gas velocity		
			Info<< " 	Filtering Up" << endl;
			sprintf(charfilterwidth, "tildeUp_%dX%dX%d",Filter,Filter,Filter);	

			volVectorField tildeUp
			(
			    IOobject
			    (
		        	charfilterwidth,
		        	runTime.timeName(),
		        	mesh,
		        	IOobject::NO_READ
			    ),
			    Up_
			);

			forAll(mesh.cells(),cellI)
			{
		    		scalar total_volume = 0;
				tildeUp[cellI] = vector(0,0,0); 			
				const labelList fcellI = stencillistfilter[filterwidth][cellI]; 

				forAll(fcellI,filtercellI)
				{
					total_volume   +=       mesh.V()[fcellI[filtercellI]];
					tildeUp[cellI] +=   (1.-voidfraction_[fcellI[filtercellI]]) 
					        	  * Up_[fcellI[filtercellI]] 
							  * mesh.V()[fcellI[filtercellI]]; 
				}
	        		if( total_volume > 0 )
				{
					tildeUp[cellI] /= total_volume; 

					if ( barvoidfraction[cellI] > 0 )
					{
						tildeUp[cellI] /= (1.-barvoidfraction[cellI]);
					}
				}
				else 
				{
					tildeUp[cellI] = vector(0,0,0);
				}
			} 

			Info<< " 	Writing filtered Up" << endl;
	    		tildeUp.write();		
			tildeUp_ = tildeUp;		
		}

    }
	
	void CalculateDragForce
	(
		cfdemCloud& sm,
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volScalarField& rho_,
		const bool& verbose_,
		vectorField& DragForce_
	)
	{
		
		// get viscosity field
		#ifdef comp
		    const volScalarField nufField = sm.turbulence().mu()/rho_;
		#else
		    const volScalarField& nufField = sm.turbulence().nu();
		#endif

		// Local variables	
		label  cellI=0;
		vector drag(0,0,0);
		scalar Volp(0);
		vector Ufluid(0,0,0);
		
		vector position(0,0,0);
		scalar voidfraction(1);
		
		vector Up(0,0,0);
		vector Ur(0,0,0);
		scalar ds(0);
		
		scalar nuf(0);
		scalar rho(0);
		scalar magUr(0);
		scalar Rep(0);
		scalar CD(0);
		scalar alps(0);
		scalar betaP(0);

		interpolationCellPoint<scalar> voidfractionInterpolator_(alpf_);
		interpolationCellPoint<vector> UInterpolator_(Uf_);	
		
		// 
		DragForce_.resize(sm.numberOfParticles());
			
		for(int index = 0; index <  sm.numberOfParticles(); index++)
	    	{
			cellI = sm.cellIDs()[index][0];
			drag = vector(0,0,0);
			Volp = 0;
			Ufluid = vector(0,0,0);
			DragForce_[index] = vector(0,0,0);
			    
			if (cellI > -1) // particle Found
			{
				position = sm.position(index);
				voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
				Ufluid = UInterpolator_.interpolate(position,cellI);

				Up = sm.velocity(index);
				Ur = Ufluid-Up;
				ds = 2*sm.radius(index);		

				nuf = nufField[cellI];
				rho = rho_[cellI];
				magUr = mag(Ur);
				Volp = ds*ds*ds*M_PI/6;
				alps = 1-voidfraction+SMALL;

				if (magUr > 0)
				{
					// calc particle Re number
					Rep = voidfraction*ds*magUr/(nuf+SMALL);

					// calc CD
					if (Rep < 1000)
					{
						CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-2.65);
					}
					else
					{
						CD = 0.44*pow(voidfraction,-2.65); 
					}  

					// calc drag coefficient 
					betaP = 3./4.*rho*CD*magUr/ds;  	

					// calc particle's drag
					drag = Volp*betaP*Ur;

				}

			}
			
			if(verbose_ && index <= 1)
			{
				Info << "" << endl;
				Info << " index = " << index << endl;
				Info << " position = " << position << endl; 
				Info << " Up = " << Up << endl;
				Info << " Ur = " << Ur << endl;
				Info << " dp = " << ds << endl;
				Info << " rho = " << rho << endl;
				Info << " nuf = " << nuf << endl;
				Info << " voidfraction = " << voidfraction << endl;
				Info << " Rep = " << Rep << endl;
				Info << " Volp = " << Volp << endl;
				Info << " alps = " << alps << endl;
				Info << " CD = " << CD << endl;
				Info << " betaP = " << betaP << endl;		   		    
				Info << " drag = " << drag << endl;
				Info << " " << endl;
			}
			
			for(int j=0;j<3;j++) DragForce_[index][j] = drag[j];		
				
		}	
	}
	
	
	
	void filteredDragCoefficient
	(
		cfdemCloud& sm,
		const bool& verbose_,
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volVectorField& tildeUf_,
		const volScalarField& barPg_,
		const volScalarField& p_,
		const volScalarField& rho_,
		const vectorField& DragForce_,
		vectorField& ResolvedDragForce_,
		vectorField& gii_
	)
	{
		
		// Local variables	
		label  cellI=0;
		vector gradPg_int(0,0,0);
		vector gradbarPg_int(0,0,0);		
		
		vector position(0,0,0);	
	    	vector Ufluid(0,0,0);
		vector tildeUfluid(0,0,0);			
		vector Up(0,0,0);
		vector Ur(0,0,0);
		
		scalar ds(0);
		scalar Volp(0);
		
		volVectorField gradp_ = fvc::grad(p_);
		volVectorField gradbarPg_ = fvc::grad(barPg_);
			   
		interpolationCellPoint<vector> gradPgInterpolator_(gradp_);
		interpolationCellPoint<vector> gradbarPgInterpolator_(gradbarPg_);
		interpolationCellPoint<vector> UInterpolator_(Uf_);	
		interpolationCellPoint<vector> tildeUInterpolator_(tildeUf_);	
		   
	    // Filtered drag coefficient 
		vectorField barBetai_(sm.numberOfParticles());
		gii_.resize(sm.numberOfParticles());
		
	     // Calculate resolved drag force
	     	CalculateDragForce
		(
			sm,
			alpf_,
			tildeUf_,
			rho_,
			verbose_,
			ResolvedDragForce_
		);  	
	
		for(int index = 0; index <  sm.numberOfParticles(); index++)
	    	{
			cellI = sm.cellIDs()[index][0];						
			if( cellI > -1 )
			{
					position = sm.position(index);
					gradPg_int = gradPgInterpolator_.interpolate(position,cellI);
					gradbarPg_int = gradbarPgInterpolator_.interpolate(position,cellI);
					Ufluid = UInterpolator_.interpolate(position,cellI);					
					tildeUfluid = tildeUInterpolator_.interpolate(position,cellI);
					Up = sm.velocity(index);
					Ur = tildeUfluid-Up;
					ds = 2*sm.radius(index);
					Volp = ds*ds*ds*M_PI/6;
					
					//Info << index << " " << gradPg_int << " " << gradbarPg_int << " " << Ufluid << " " << " " << tildeUfluid << endl;
 										
					for(int j=0;j<3;j++) 
					{						
							barBetai_[index][j]  = ( - Volp * gradPg_int[j]
								 		 + Volp * gradbarPg_int[j]
										 + DragForce_[index][j] 
										) / ( Ur[j] + SMALL ) ;
							
							gii_[index][j] = barBetai_[index][j] * (Ufluid[j] - Up[j])
								            / ( DragForce_[index][j] + SMALL );
									    
							if(verbose_ && index <= 1)
							{
								Info << " barBetai_[" << index <<"]["<<j<<"] = " << barBetai_[index][j] << endl; 
								Info << "      gii_[" << index <<"]["<<j<<"] = " << gii_[index][j]      << endl; 
							}
							
							
					}	
			}	

				
		}
		
		//Info << gii_ << endl ;		
		 
   	}	
	 
 	void createParcels
	(
		const fvMesh& mesh,			
		cfdemCloud& sm,
		cfdemCloud& parcelsm,
		const int& nparticle_,
		const bool verbose_,
		const fileName outputRelativePath_, // Write parcel diameter in a file in the directory of ./PostProcessing..	
		labelListList& parcelList,
		const bool& weighting_	 			
	)
	{		

		if ( nparticle_ > sm.numberOfParticles() )
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
		nPts =  sm.numberOfParticles();		
		
		Info << " Number of particles in a parcel = " << nparticle_ << endl;
		
		// Squared radius distance between nearest neighbours: defined as max. distance in the domain
		const pointField& pp = mesh.points();

		// Min, max x-coordinates	
		scalar min_x = Foam::min(pp & vector(1,0,0));
		scalar max_x = Foam::max(pp & vector(1,0,0));

		// Min, max y-coordinates		
		scalar min_y = Foam::min(pp & vector(0,1,0));
		scalar max_y = Foam::max(pp & vector(0,1,0));

		// Min, max z-coordinates		
		scalar min_z = Foam::min(pp & vector(0,0,1));
		scalar max_z = Foam::max(pp & vector(0,0,1));		
		
		/*	
		const scalar sqRad =  ( (max_x - min_x) * (max_x - min_x) )
				     +( (max_y - min_y) * (max_y - min_y) )
				     +( (max_z - min_z) * (max_z - min_z) );	
		*/
                      
                const scalar sqRad = ( ( (max_x - min_x) * (max_x - min_x) )
                                      +( (max_y - min_y) * (max_y - min_y) )
                                      +( (max_z - min_z) * (max_z - min_z) ) ) / 10.;    
		
		Info << " Squared radius = " << sqRad << endl;
		
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
		labelList particleList(sm.numberOfParticles());				
		for(int index = 0; index <  sm.numberOfParticles(); index++)
	        {			
				dataPts[index][0] = sm.position(index).x();		
				dataPts[index][1] = sm.position(index).y();		
				dataPts[index][2] = sm.position(index).z();
				particleList[index] = index;			
		}
				
		kdTree = new ANNkd_tree(dataPts, nPts, dim);
		
		label parcelI = 0;						
		forAll(particleList,index)
		{				
						
			if ( particleList[index] > -1 )
				{					
					
					//if (verbose_) Info << " Parcel = " << parcelI << endl;
					
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
							
					//if (verbose_) Info << "NN:  Index  Distance\n";
										
					int sumnnIdx = 0; 
 	
					for (int i = 0; i < k; i++)
					{
						if ( particleList[nnIdx[i]] == -1 )
						{
							
							particleList[nnIdx[i]] = sm.numberOfParticles() + 1;
								         sumnnIdx += particleList[nnIdx[i]];
						}
						else
						{
 							
							particleList[nnIdx[i]] = -1 ;
								         sumnnIdx += particleList[nnIdx[i]];
							parcelList[parcelI][i] = nnIdx[i] ;
							
						}
						
						//if (verbose_ ) Info << index << " " << i << " " << nnIdx[i] << " " << particleList[nnIdx[i]]<< " " << dists[i] << "\n";
															 						 
					}
										
					if ( sumnnIdx == -nparticle_ )
					{
						parcelI++;																				
					}
					
				}	
				
		}
				
		Info << " Number of parcels = " << parcelI << endl;
		Info << " Number of particles not used to construct parcels = " << sm.numberOfParticles()-parcelI*nparticle_ << endl;
				
		// Resize parcel list 
		parcelList.resize(parcelI);	
		
		// Parcel center of mass
		vectorField parcelCenter(parcelI,vector(0,0,0));
		// Calculate weighted velocities
		vectorField parcelVel(parcelI,vector(0,0,0));
		// Parcel diameter
		scalarField parcelDiameter(parcelI,scalar(0.0));	
		
		// Set number of parcels
		parcelsm.setNumberOfParticles(parcelI);	
						
		vector Up(0,0,0);							
		for(int parcelII = 0; parcelII < parcelsm.numberOfParticles(); parcelII++)
		{ 			
			labelList parcelListL = parcelList[parcelII];
			Up = vector(0,0,0);
			parcelDiameter[parcelII] = 0;
						
			forAll(parcelListL,particleI)
			{				
				
		
				Up = sm.velocity(parcelListL[particleI]);
									
				for(int j=0;j<3;j++) 
				{	
						// Parcel velocity 
						parcelVel[parcelII][j]    += Up[j];
						// Parcel center of mass
						parcelCenter[parcelII][j] += sm.position(parcelListL[particleI])[j];	
		        	      					
				}	
				parcelDiameter[parcelII] =  2.*nparticle_*sm.radius(parcelListL[particleI]);											
			}
			
			//Info << " Parcel I " << parcelII << " Parcel center " <<  parcelCenter[parcelII] << endl;
				
		}			    
				
		// Define Parcel Cloud positions, velocities, radius and cellIDs				
		double **parcelPositions;
		double **parcelVelocities;
		double **parcelRadii;
		double **parcelCellID;

		parcelsm.dataExchangeM().allocateArray(parcelPositions,0.,3);
		parcelsm.dataExchangeM().allocateArray(parcelVelocities,0.,3);
		parcelsm.get_radii(parcelRadii); 
		parcelsm.get_cellIDs(parcelCellID);				
		
		for(int index = 0; index <  parcelsm.numberOfParticles(); index++)
	    	{
			for(int j=0;j<3;j++) 
			{
				parcelPositions[index][j]  = parcelCenter[index][j]/nparticle_;
				parcelVelocities[index][j] = parcelVel[index][j];
			}
			parcelRadii[index][0] = parcelDiameter[index]/2.; 
			
		}
		
		parcelsm.locateM().findCell(NULL,parcelPositions,parcelCellID,parcelsm.numberOfParticles());
		parcelsm.setPos(parcelPositions);
		parcelsm.setVel(parcelVelocities);	
				
		if ( verbose_ && nparticle_ <=2)
		{							
			int index = 0;
			labelList parcelPart = parcelList[index];
			Info << " Parcel particle list " << parcelPart << endl; 
			Info << " Parcel center     " <<  parcelsm.position(index)  		<< endl;	
			Info << " Parcel velocity   " <<  parcelsm.velocity(index)  		<< endl;
			Info << " Parcel diameter   " << 2.*parcelsm.radius(index) 			<< endl;
					
			forAll(parcelPart, ii)
			{
				Info << " Particle " << parcelPart[ii] << endl;
				Info << " Particle center    " <<        sm.position(parcelPart[ii])   << endl;
				Info << " Particle velocity  " <<        sm.velocity(parcelPart[ii])   << endl;
				Info << " Particle diameter  " <<       2.*sm.radius(parcelPart[ii]) 	<< endl;
			}			
		}
										
	}
 

 	void ParticleCoarsening
	(
		const fvMesh& mesh,
		const labelListListList& stencillistfilter, // Used to check the parcel diameter > filter width
		const int& filterwidth,				
		cfdemCloud& sm,
		cfdemCloud& parcelsm,
		const bool verbose_,
		const volVectorField& Uf_,
		const volVectorField& tildeUf_,
		const volScalarField& p_,
		const volScalarField& barPg_,
		const vectorField& DragForce_,
		const vectorField& ResolvedDragForce_,
		vectorField& parcelgii_,
		const fileName outputRelativePath_, // Write parcel diameter in a file in the directory of ./PostProcessing..		
		const labelListList& parcelList,
		const bool& weighting_
	)
 	{
 		// Distance vector
		vector dist(0,0,0);
		
		// Parcel averaged force 
		vectorField parcelForce(parcelsm.numberOfParticles(),vector(0,0,0));
		
		// Parcel averaged drag force
		vectorField parcelDragForce(parcelsm.numberOfParticles(),vector(0,0,0));
		
		// Parcel averaged filtered relative velocity 
		vectorField parcelFilteredRelativeVelocity(parcelsm.numberOfParticles(),vector(0,0,0));
		
		// Parcel averaged relative velocity
		vectorField parcelRelativeVelocity(parcelsm.numberOfParticles(),vector(0,0,0));		
		
		volVectorField gradp_ = fvc::grad(p_);
		volVectorField gradbarPg_ = fvc::grad(barPg_);
			   
		interpolationCellPoint<vector> gradPgInterpolator_(gradp_);
		interpolationCellPoint<vector> gradbarPgInterpolator_(gradbarPg_);
		
		interpolationCellPoint<vector> UInterpolator_(Uf_);	
		interpolationCellPoint<vector> tildeUInterpolator_(tildeUf_);
		
		vector gradPg_int(0,0,0);	
	   	vector gradbarPg_int(0,0,0);		
		vector position(0,0,0);
		vector Up(0,0,0);
		vector Ufluid(0,0,0);
		vector tildeUfluid(0,0,0);
		vector Ur(0,0,0);
		scalar Volp(0);
		scalar ds(0);
		
		label cellI(0);
					
		//forAll(parcelList,parcelII)	
		for(int parcelII = 0; parcelII <  parcelsm.numberOfParticles(); parcelII++)
		{ 			
			labelList parcelListL = parcelList[parcelII];
			dist = vector(0,0,0);
			
			//Info << " Parcel   " << parcelII << endl;
			
			forAll(parcelListL,particleI)
			{				
				//Info << " Particle " << parcelListL[particleI] << endl;
				
				if(weighting_)
				{
					dist[0] = sqrt(   ( sm.position(parcelListL[particleI]).x() - parcelsm.position(parcelII).x() ) 
					        	* ( sm.position(parcelListL[particleI]).x() - parcelsm.position(parcelII).x() ) );

					dist[1] = sqrt(   ( sm.position(parcelListL[particleI]).y() - parcelsm.position(parcelII).y() ) 
					        	* ( sm.position(parcelListL[particleI]).y() - parcelsm.position(parcelII).y() ) );

					dist[2] = sqrt(   ( sm.position(parcelListL[particleI]).z() - parcelsm.position(parcelII).z() ) 
					        	* ( sm.position(parcelListL[particleI]).z() - parcelsm.position(parcelII).z() ) );				
				}
				
				cellI = sm.cellIDs()[parcelListL[particleI]][0];
				gradPg_int = vector(0,0,0);
				gradbarPg_int = vector(0,0,0);
				Ufluid = vector(0,0,0);
				Up = vector(0,0,0);
				Ur = vector(0,0,0);
				Volp = 0;
				tildeUfluid = vector(0,0,0);
											
				if( cellI > -1 )
				{
					position = sm.position(parcelListL[particleI]); 
					gradPg_int = gradPgInterpolator_.interpolate(position,cellI); 
					gradbarPg_int = gradbarPgInterpolator_.interpolate(position,cellI); 
					Ufluid = UInterpolator_.interpolate(position,cellI); 
					tildeUfluid = tildeUInterpolator_.interpolate(position,cellI);	
					Up = sm.velocity(parcelListL[particleI]);
					Ur = tildeUfluid-Up;					 	
				    	ds = 2*sm.radius(parcelListL[particleI]);
					Volp = ds*ds*ds*M_PI/6;
																				
					// Parcel averaged resolved force
					for(int j=0;j<3;j++) 
					{		

						if(weighting_)
						{
							            parcelForce[parcelII][j] += ( - Volp * gradPg_int[j] 
							                                          + Volp * gradbarPg_int[j]
								 		 	          + DragForce_[parcelListL[particleI]][j]   ) / ( dist[j] + SMALL ) ; 
										
						 parcelFilteredRelativeVelocity[parcelII][j] += Ur[j] / ( dist[j] + SMALL); 
					
						                parcelDragForce[parcelII][j] += DragForce_[parcelListL[particleI]][j] / ( dist[j] + SMALL );
							
						         parcelRelativeVelocity[parcelII][j] += (Ufluid[j] - Up[j]) / ( dist[j] + SMALL ); 
						}
						
						else

						{
						
							            parcelForce[parcelII][j] += ( - Volp * gradPg_int[j] 
							                                          + Volp * gradbarPg_int[j]
								 		 	          + DragForce_[parcelListL[particleI]][j]   )  ; 
										
						 parcelFilteredRelativeVelocity[parcelII][j] += Ur[j] ; 
					
						                parcelDragForce[parcelII][j] += DragForce_[parcelListL[particleI]][j] ;
							
						         parcelRelativeVelocity[parcelII][j] += (Ufluid[j] - Up[j]) ; 						
						}
							 										
						
					}					

				}
				
					
			}										
			
		}				

		// Parcel averaged filtered drag coefficient
		vectorField parcelbarBetai_(parcelsm.numberOfParticles(),vector(0,0,0));
		
		// Drag correction with particle coarsening
		parcelgii_.resize(parcelsm.numberOfParticles());
		
		forAll(parcelList,parcelII)
		{
			for(int j=0;j<3;j++) 
			{				
			 parcelbarBetai_[parcelII][j] =        parcelForce[parcelII][j] / ( parcelFilteredRelativeVelocity[parcelII][j] + SMALL ) ;
			      parcelgii_[parcelII][j] =    parcelbarBetai_[parcelII][j] *           parcelRelativeVelocity[parcelII][j] 
				                       / ( parcelDragForce[parcelII][j] + SMALL ) ;				
			}
		}		
		
		if( verbose_ )
		{
			int index = 0;
			Info << " Parcel barBetai " << parcelbarBetai_[index] << endl;
			Info << " Parcel gii      " << parcelgii_[index]      << endl;		
		}
	
	}
	
	void applyBins
	(
		cfdemCloud& sm,
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
		const volScalarField& baralpf_,
		const volVectorField& tildeUf_,
		vectorField& gii_,
		Field< Field < Field <scalar> > >& 	giiCondalppFilter_,
		Field< Field < Field <scalar> > >& 	velSlipCondalppFilter_,
		Field < Field <scalar> >& 		taupVelSlipCondalppFilter_,
		Field < Field <scalar> >& 		numberCondalppFilter_, 		
		Field< Field < Field <scalar> > >&	VelSlipJointgiiCondalppFilter_,
		Field< Field < Field <scalar> > >&	TaupVelSlipJointgiiCondalppFilter_,	
		Field <Field < Field <scalar> > >& 	numberVelSlipJointgiiCondalppFilter_,
		Field <Field < Field <scalar> > >& 	numberTaupVelSlipJointgiiCondalppFilter_
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
			FatalError << " Max. solid volumde fraction equal to = " << maxalpp_ << abort(FatalError);
		}	
		
		// Binning particle solid volume fraction
		const scalar binDx_   = (maxalpp_-minalpp_)/nBin_;	
		// Binning relative velocity 
		const scalar binVelSlip_ = (maxVelSlip_-minVelSlip_)/nBin_;
		// Binning tau_p |Vr|
		const scalar binTaupVelSlip_= (maxTaupVelSlip_-minTaupVelSlip_)/nBin_;
		
		interpolationCellPoint<scalar> baralpfInterpolator_(baralpf_);		
		interpolationCellPoint<vector> tildeUInterpolator_(tildeUf_);
		
		label cellI = 0;
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
		
		for(int index = 0; index <  sm.numberOfParticles(); index++)
	    	{
			cellI = sm.cellIDs()[index][0];
									    
			if (cellI > -1) // particle Found
			{
								
				Up            = sm.velocity(index);
				position      = sm.position(index);
				baralpf       =    baralpfInterpolator_.interpolate(position,cellI) ;
				baralpp       = 1.-baralpfInterpolator_.interpolate(position,cellI) ;
				intertildeUf_ =     tildeUInterpolator_.interpolate(position,cellI) ;
			
				Ur = intertildeUf_-Up;
				ds = 2*sm.radius(index);		
				nuf = nufField[cellI];
				rhof = rho_[cellI];
				magUr = mag(Ur);
								
				label binI   = floor((baralpp-minalpp_)/binDx_);
								
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
			
                                if ( binI < nBin_) // Avoid to sample out of the defined binning interval
				{					
					for(int j=0;j<3;j++) 
					{	
					      giiCondalppFilter_[binI][j][filterwidth] += gii_[index][j];
				          velSlipCondalppFilter_[binI][j][filterwidth] += intertildeUf_[j]-Up[j];  										
					}
				         taupVelSlipCondalppFilter_[binI][filterwidth] += taup;				      					      					      
					      numberCondalppFilter_[binI][filterwidth] += 1;				
								
				}
				
				// Along z-direction
				label binVelSlipI     = floor((intertildeUf_[2]-Up[2]-minVelSlip_)/binVelSlip_); 
																
				label binTaupVelSlipI = floor((taup*magUr-minTaupVelSlip_)/binTaupVelSlip_);
							
				if ( binI < nBin_ && binVelSlipI < nBin_ && binVelSlipI > 0 )
				{
				   	         VelSlipJointgiiCondalppFilter_[binI][binVelSlipI][filterwidth] += gii_[index][2];
					   numberVelSlipJointgiiCondalppFilter_[binI][binVelSlipI][filterwidth] += 1;
				} 
				
				if ( binI < nBin_ && binTaupVelSlipI < nBin_ )
				{    
					 TaupVelSlipJointgiiCondalppFilter_[binI][binTaupVelSlipI][filterwidth] += gii_[index][2];
				   numberTaupVelSlipJointgiiCondalppFilter_[binI][binTaupVelSlipI][filterwidth] += 1;
				}
				
					
			}
		}		
	} 

	void writeBins
	(
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
		bool verbose_	
	)		
	{

		// Create output folder	
		if( !isDir(mesh.time().path()/outputRelativePath_) )
		{
			mkDir(mesh.time().path()/outputRelativePath_);
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
			
			if ( nparticle_ > 0 )
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
			
			OFstream         str_gii(mesh.time().path()/outputRelativePath_/outputfile);
			OFstream     str_velSlip(mesh.time().path()/outputRelativePath_/outputfileVelSlip);
			OFstream str_TaupVelSlip(mesh.time().path()/outputRelativePath_/outputfileTaupVelSlip);			
			OFstream       str_nreal(mesh.time().path()/outputRelativePath_/outputfileN);

			OFstream                str_JointVelSlip(mesh.time().path()/outputRelativePath_/outputfileJointVelSlip);
			OFstream            str_JointTaupVelSlip(mesh.time().path()/outputRelativePath_/outputfileJointTaupVelSlip);			
			OFstream           str_JointNrealVelSlip(mesh.time().path()/outputRelativePath_/outputfileJointNVelSlip);
			OFstream       str_JointNrealTaupVelSlip(mesh.time().path()/outputRelativePath_/outputfileJointNTaupVelSlip);

			
			Info << " 		   " << endl; 		
			Info << " Writing the file " << outputfile << ", " << outputfileVelSlip      << ", " << outputfileTaupVelSlip      << ", " << outputfileN 
						     		   << ", " << outputfileJointVelSlip << ", " << outputfileJointTaupVelSlip << ", " << outputfileJointNVelSlip << ", " <<  outputfileJointNTaupVelSlip
			     << " into the folder: " << mesh.time().path()/outputRelativePath_ << endl;
			
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
					Info  	<< minalpp_+(i+1./2.)*(maxalpp_-minalpp_)/nBin_   
					 << " " << giiCondalppFilter_[i][0][filterwidth] 
					 << " " << giiCondalppFilter_[i][1][filterwidth]  
					 << " " << giiCondalppFilter_[i][2][filterwidth] 
					 << " " << numberCondalppFilter_[i][filterwidth] 						<< endl;
	
					Info  	<< (minalpp_+i+1./2.)*(maxalpp_-minalpp_)/nBin_   
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	

	Foam::timeSelector::addOptions();
	#   include "addRegionOption.H"
	Foam::argList::addBoolOption
	(
	"noWrite",
	"suppress writing results"
	);
	#include "addDictOption.H"

	#include "setRootCase.H"
	#include "createTime.H"
	Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
	#include "createNamedMesh.H"

	#include "createFields.H" 	 

	// Create particle cloud			
	cfdemCloud particleCloud(mesh);
	
	particleCloud.reAllocArrays();

	double **positions;
	double **velocities;
	double **radii;
	double **cellID;

	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
	
	
	// Create parcel cloud	( Max:5 )
	cfdemCloud parcelCloud(mesh);
	parcelCloud.reAllocArrays();
	
	cfdemCloud parcelCloud_2(mesh);
	parcelCloud_2.reAllocArrays();
	
	cfdemCloud parcelCloud_3(mesh);
	parcelCloud_3.reAllocArrays();
	
	cfdemCloud parcelCloud_4(mesh);
	parcelCloud_4.reAllocArrays();
	
	cfdemCloud parcelCloud_5(mesh);
	parcelCloud_5.reAllocArrays();					
		
	// Parcel velocities are defined in parcel coarsening subroutine	
	
	// Post-processing filtered drag dictionary 
	// Read dictionary + dictionaryProps 
	const dictionary dict(particleCloud.couplingProperties());
	// Data exchange oneWayVTK
	/*const dictionary VTKDict(dict.subDict("oneWayVTKProps"));					
	const scalar DEMts(readScalar(VTKDict.lookup("DEMts")));
    	fileName relativePath(VTKDict.lookup("relativePath"));
	fileName couplingFilename(VTKDict.lookup("couplingFilename"));
	const int maxNumberOfParticles(readScalar(VTKDict.lookup("maxNumberOfParticles")));	*/
				
	// Wen&Yu Drag
	const dictionary propsDict(dict.subDict("PostProcessingFilteredDrag"));					
	bool verbose(false);
	//const int DEM_dump_Interval(readScalar(propsDict.lookup("DEM_dump_Interval")));	
	const int nBin(readScalar(propsDict.lookup("nBin")));
	const scalar maxalpp(readScalar(propsDict.lookup("Max_alpp")));
	const scalar minalpp(readScalar(propsDict.lookup("Min_alpp")));
	// Binning with VelSlip, taup |Vr|
	//const scalar minVelSlip(-1);
	//const scalar maxVelSlip( 1);
	//const scalar minTaupVelSlip(0);
	//const scalar maxTaupVelSlip(0.5);
	//
	const scalar minVelSlip(readScalar(propsDict.lookup("Min_VelSlip")));
	const scalar maxVelSlip(readScalar(propsDict.lookup("Max_VelSlip")));
	const scalar minTaupVelSlip(readScalar(propsDict.lookup("Min_TaupVelSlip")));
	const scalar maxTaupVelSlip(readScalar(propsDict.lookup("Max_TaupVelSlip")));	
	//
	const scalar rhop(readScalar(propsDict.lookup("rhopart")));
    	fileName outputRelativePath(propsDict.lookup("outputRelativePath"));
	fileName stencilsOutputRelativePath(propsDict.lookup("stencilsOutputRelativePath"));
	bool particleCoarsening(false);
	bool verboseParticleCoarsening(false);
	bool weighting(false);
	//int npartParticleCoarsening(1);
	labelList npartParticleCoarsening(1);
	if(propsDict.found("verbose")) verbose = true;	
	if(propsDict.found("ParticleCoarsening")) 
	{
		particleCoarsening = true;
		if(propsDict.found("verboseParticleCoarsening")) verboseParticleCoarsening = true; 
		npartParticleCoarsening = propsDict.lookup("npartParticleCoarsening");
		labelList test(propsDict.lookup("npartParticleCoarsening"));
		if ( test.size() > 5 )  FatalError << " Parcel list size cannot be greater than 5 " << abort(FatalError);
	}
	
	if(propsDict.found("weighting")) weighting = true;
	
	// Filtering Euler-Euler Approach
	bool EulerianParticleVelocity(false);
	if(propsDict.found("EulerianParticleVelocity")) EulerianParticleVelocity = true;
	
	// Filter parameters
	#include "FilterVariables.H"	

	Foam::constructfilter //_ver2
	(
		args,
		runTime,
		mesh,
		minFilterWidth,
		maxFilterWidth,
		FilterIncrement,
		StencilListFilter,
		stencilsOutputRelativePath
	);				
		
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();

		if ( runTime.timeName() != "0" )
		{
			Foam::readEulerianVariables
			(
				args, 
				runTime, 
				mesh,
				voidfraction,
				U,
				rho,
				p
			);
			
			int count = runTime.value() / particleCloud.dataExchangeM().DEMts();				
			// timeI > 0 (particles are not initialised in the folder "0")

			Info<< " " << endl;
			particleCloud.dataExchangeM().getData("v","vector-atom",velocities,count);
			Info<< " 	Reading particle velocities" << endl;
			Info<< " " << endl;
			
			particleCloud.dataExchangeM().getData("x","vector-atom",positions,count);
			Info<< " 	Reading particle positions" << endl;
			Info<< " " << endl;

			particleCloud.dataExchangeM().getData("radius","scalar-atom",radii,count);
			Info<< " 	Reading particle radius" << endl;		
			Info<< " " << endl;
		
			particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
			particleCloud.setPos(positions);
			particleCloud.setVel(velocities);

			if(verbose)
			{
				int index = 0;
				Info << "" << endl;
				Info << " index  = " << index << endl;
				Info << " rp     = " << particleCloud.radius(index) << endl;
				Info << " Vp     = " << particleCloud.velocity(index) << endl;
				Info << " Xp     = " << particleCloud.position(index) << endl;
				Info << " CellID = " << particleCloud.particleCell(index) << endl;
				Info << " Ug     = " << U[particleCloud.particleCell(index)] << endl;
				Info << "" << endl;
			}

			// 	Call DragForce function for the drag force	
			Foam::CalculateDragForce
			(
				particleCloud,
				voidfraction,
				U,
				rho,
				verbose,
				DragForce
			);

			if(EulerianParticleVelocity)
			{
				
				
				Foam::EulerianParticleVelocity
				(
					particleCloud,			
					mesh,
					Up			
				);
								
				Up.write();								
				
				volVectorField EulerianDragForce = U-Up;				
				Foam::CalculateEulerianDragForce
				(
					particleCloud,
					mesh,
					voidfraction,
					U,
					Up,
					rho,
					EulerianDragForce	
				);

			}


			// Parcel particle list
			labelListListList parcelParticleListList;
			
			if(particleCoarsening)
			{
				//forAll(npartParticleCoarsening,npartI)
				for ( int npartI = 0; npartI < npartParticleCoarsening.size(); npartI++ )				
				{
					Info << " " << endl;
					Info << " Particle coarsening: "<< endl;
					
					// Number of parcel
					int nParcel =  particleCloud.numberOfParticles()/npartParticleCoarsening[npartI];					
					labelList particleID(npartParticleCoarsening[npartI],(0));
					labelListList parcelParticleList(nParcel,particleID);
					//parcelParticleListList.setSize(5,parcelParticleList);		
					parcelParticleListList.setSize(npartParticleCoarsening.size(),parcelParticleList);
					
					//Info << " parcelParticleListList " << parcelParticleListList[npartI] << endl;					
					
					if ( npartI == 0 )
					Foam::createParcels(mesh,particleCloud,parcelCloud,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										
					
										parcelParticleListList[npartI] = parcelParticleList;
										
					if ( npartI == 1 )
					Foam::createParcels(mesh,particleCloud,parcelCloud_2,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										

										parcelParticleListList[npartI] = parcelParticleList;
					
					if ( npartI == 2 )
					Foam::createParcels(mesh,particleCloud,parcelCloud_3,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										

										parcelParticleListList[npartI] = parcelParticleList;

					if ( npartI == 3 )
					Foam::createParcels(mesh,particleCloud,parcelCloud_4,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										

										parcelParticleListList[npartI] = parcelParticleList;
										
					if ( npartI == 4 )
					Foam::createParcels(mesh,particleCloud,parcelCloud_5,npartParticleCoarsening[npartI],
										verboseParticleCoarsening,outputRelativePath,parcelParticleList,weighting);										

										parcelParticleListList[npartI] = parcelParticleList;		
																		
					//Info << " parcelParticleListList " << parcelParticleListList[npartI] << endl;	
										
				}
			}							
															
			for(int FilterWidth = minFilterWidth; FilterWidth <= maxFilterWidth; FilterWidth += FilterIncrement )
			{
				Info << " " << endl;
				int Filter = 2*(FilterWidth+1)+1;
				Info << "Filter size = " << Filter << "X" << Filter << "X" << Filter << endl;

				Foam::filteringEulerianVariables
				(
					args, 
					runTime, 
					mesh,
					StencilListFilter,
					FilterWidth,
					voidfraction,
					U,
					p,
					baralpf,				
					tildeUf,
					barPg,
					Up,
					tildeUp,
					EulerianParticleVelocity
				);							

				Foam::filteredDragCoefficient
				(
					particleCloud,
					verbose,
					baralpf,
					U,
					tildeUf,
					barPg,
					p,
					rho,
					DragForce,
					ResolvedDragForce,
					gii
				);

				if(particleCoarsening)
				{
					//forAll(npartParticleCoarsening,npartI)
					for ( int npartI = 0; npartI < npartParticleCoarsening.size(); npartI++ )
					{
						// Particle coarsening active or not
						if ( npartI == 0 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
																
							Foam::applyBins(parcelCloud,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );
									

																
						}	
						
						if ( npartI == 1 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud_2,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
								
							Foam::applyBins(parcelCloud_2,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );						
									
						}
						
						if ( npartI == 2 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud_3,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
								
							Foam::applyBins(parcelCloud_3,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );

						}
						
						if ( npartI == 3 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud_4,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
								
							Foam::applyBins(parcelCloud_4,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );

						}
						
						if ( npartI == 4 )
						{
							Foam::ParticleCoarsening(mesh,StencilListFilter,FilterWidth,particleCloud,parcelCloud_5,
								verboseParticleCoarsening,U,tildeUf,p,barPg,DragForce,ResolvedDragForce,
								parcelgii,outputRelativePath,parcelParticleListList[npartI],weighting);
								
							Foam::applyBins(parcelCloud_5,
									FilterWidth,maxFilterWidth,
									nBin,
									minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
									rhop,rho,
									baralpf,
									tildeUf,
									parcelgii,
									NparcelgiiCondalppFilter[npartI],
									NparcelVelSlipCondalppFilter[npartI],
									NparcelTaupVelSlipCondalppFilter[npartI],
									NparcelnumberCondalppFilter[npartI],
									NparcelVelSlipJointgiiCondalppFilter[npartI],
									NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
									NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI]       );								

								     
						}
						
					}
				}							

				Foam::applyBins
				(				
					particleCloud,
					FilterWidth,
					maxFilterWidth,
					nBin,
					minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,
					rhop,rho,
					baralpf,
					tildeUf,
					gii,
					giiCondalppFilter,					
					velSlipCondalppFilter,
					taupVelSlipCondalppFilter,
					numberCondalppFilter,
					VelSlipJointgiiCondalppFilter,
					TaupVelSlipJointgiiCondalppFilter,
					numberVelSlipJointgiiCondalppFilter,
					numberTaupVelSlipJointgiiCondalppFilter				
				);						
			
			}	
		
		}	

	}
	
	
	// Write bins		
	Foam::writeBins
	(
		mesh,
		nBin,
		minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,		
		outputRelativePath,
		minFilterWidth, 
		maxFilterWidth,
		FilterIncrement,
		-1,				// Just for fluid coarsening results
		giiCondalppFilter,
		velSlipCondalppFilter,
		taupVelSlipCondalppFilter,
		numberCondalppFilter,
		VelSlipJointgiiCondalppFilter,
		TaupVelSlipJointgiiCondalppFilter,
		numberVelSlipJointgiiCondalppFilter,
		numberTaupVelSlipJointgiiCondalppFilter,
		verbose						
	);
	
	if(particleCoarsening)
	{
		//forAll(npartParticleCoarsening,npartI)
		for ( int npartI = 0; npartI < npartParticleCoarsening.size(); npartI++ )		
		{
			Foam::writeBins
			(
				mesh,
				nBin,
				minalpp,maxalpp,minVelSlip,maxVelSlip,minTaupVelSlip,maxTaupVelSlip,		
				outputRelativePath,
				minFilterWidth, 
				maxFilterWidth,
				FilterIncrement,
				npartParticleCoarsening[npartI],
				NparcelgiiCondalppFilter[npartI],
				NparcelVelSlipCondalppFilter[npartI],
				NparcelTaupVelSlipCondalppFilter[npartI],
				NparcelnumberCondalppFilter[npartI],
				NparcelVelSlipJointgiiCondalppFilter[npartI],
				NparcelTaupVelSlipJointgiiCondalppFilter[npartI],
				NparcelnumberVelSlipJointgiiCondalppFilter[npartI],
				NparcelnumberTaupVelSlipJointgiiCondalppFilter[npartI],				
				verboseParticleCoarsening						
			);			
		}
	}					

		
	particleCloud.dataExchangeM().destroy(positions,3);
	particleCloud.dataExchangeM().destroy(velocities,3);

	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
		
    return 0;
}


// ************************************************************************* //
