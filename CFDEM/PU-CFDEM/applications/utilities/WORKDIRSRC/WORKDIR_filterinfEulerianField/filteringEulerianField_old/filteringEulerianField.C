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
#include "./globalIndexStencils/CPCCellToCellStencil_mod.H"
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

		Info << " Domain " << " x=[" << min_x << ":" << max_x << "]" << " y=[" << min_y << ":" << max_y << "]" << " z=[" << min_z << ":" << max_z << "]" << endl;

		const scalar cell_volume_ref = mesh.V()[0];
		forAll(mesh.cells(),cellI)
		{
			if ( abs(cell_volume_ref-mesh.V()[cellI]) > SMALL)  FatalError<< "Non-uniform mesh !!! " << abort(FatalError);
		}

		const scalar delta = pow(cell_volume_ref,1./3);
		Info << " Uniform mesh,  volume of a cell = " << cell_volume_ref << endl;
		Info << " Delta_x = Delta_y = Delta_z = " << delta << endl;


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
		Info << " Number of coupled patches = " << nCoupledPatches << endl;
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
				
				/*
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
				*/
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
		volScalarField& barPg_ 
	)
    {

		char charfilterwidth[100]; 
		
		// Filtering volume fraction of gas phase		
		Info<< " 	Filtering voidfraction" << endl;
		sprintf(charfilterwidth, "barvoidfraction_%dX%dX%d",filterwidth,filterwidth,filterwidth);	

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
		sprintf(charfilterwidth, "tildeU_%dX%dX%d",filterwidth,filterwidth,filterwidth);	
		
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

    }
	
	void DragForce
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
	     	DragForce
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
		labelListList& parcelList	 			
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
		
		// Number of parcel
		const int nParcel =  sm.numberOfParticles()/nparticle_;
		
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
			
		const scalar sqRad =  ( (max_x - min_x) * (max_x - min_x) )
				     +( (max_y - min_y) * (max_y - min_y) )
				     +( (max_z - min_z) * (max_z - min_z) );	
		
		
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

		// Create parcel list
		labelList particleID(nparticle_,(0));
		//labelListList parcelList(nParcel,particleID);
		parcelList.resize(nParcel,particleID);
		
		// Parcel center of mass
		vectorField parcelCenter(nParcel);
		label parcelI = 0;		
		//for(int index = 0; index <  sm.numberOfParticles(); index++)
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
				        parcelCenter[parcelI][0] = 0;
					parcelCenter[parcelI][1] = 0;
					parcelCenter[parcelI][2] = 0;
 	
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
						
						// Parcel center of mass
						parcelCenter[parcelI][0] += sm.position(nnIdx[i]).x();
						parcelCenter[parcelI][1] += sm.position(nnIdx[i]).y();
						parcelCenter[parcelI][2] += sm.position(nnIdx[i]).z();
						 
						 
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
		parcelList.resize(parcelI,particleList);
		parcelCenter.resize(parcelI);
		
		// Calculate parcel center of mass
		parcelCenter /= nparticle_;	
		
		// Calculate weighted velocities
		vectorField parcelVel(parcelI);		

		// Total distance of particles to center of mass
		vectorField tot_dist(parcelI);
				
		vector dist(0,0,0);
		vector Up(0,0,0);					
		// Parcel diameter
		scalarField parcelDiameter(parcelI,scalar(0.0));
		forAll(parcelList,parcelII)
		{ 			
			labelList parcelListL = parcelList[parcelII];
			dist = vector(0,0,0);
			Up = vector(0,0,0);
			parcelDiameter[parcelII] = 0;
			forAll(parcelListL,particleI)
			{				
				dist[0] = sqrt(   ( sm.position(parcelListL[particleI]).x() - parcelCenter[parcelII].x() ) 
					            * ( sm.position(parcelListL[particleI]).x() - parcelCenter[parcelII].x() ) );
				
				dist[1] = sqrt(   ( sm.position(parcelListL[particleI]).y() - parcelCenter[parcelII].y() ) 
					            * ( sm.position(parcelListL[particleI]).y() - parcelCenter[parcelII].y() ) );
				
				dist[2] = sqrt(   ( sm.position(parcelListL[particleI]).z() - parcelCenter[parcelII].z() ) 
					            * ( sm.position(parcelListL[particleI]).z() - parcelCenter[parcelII].z() ) );
						    
				if ( mag(dist) > parcelDiameter[parcelII] )
				{
					parcelDiameter[parcelII] = 2.*mag(dist);
				}

				Up = sm.velocity(parcelListL[particleI]);
				
				for(int j=0;j<3;j++) 
				{	
				      parcelVel[parcelII][j] += Up[j] / ( dist[j] + SMALL );
		        	       tot_dist[parcelII][j] +=  1./ ( dist[j] + SMALL );				
				}
			}
				
		}				    
		
		// Define Parcel Cloud positions, velocities, radius and cellIDs		
		
		// Set number of parcels
		parcelsm.setNumberOfParticles(parcelI);	
		
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
				parcelPositions[index][j]  = parcelCenter[index][j];
				parcelVelocities[index][j] = parcelVel[index][j]/tot_dist[index][j];

				if ( parcelDiameter[index] < SMALL ) // case npart=1, avoid if one define npart=1
				{			
					parcelRadii[index][j] = nparticle_*sm.radius(index);
				}
				else
				{
					parcelRadii[index][j] = parcelDiameter[index]/2.;	
				}

			}			

		}
		
		parcelsm.locateM().findCell(NULL,parcelPositions,parcelCellID,parcelsm.numberOfParticles());
		parcelsm.setPos(parcelPositions);
		parcelsm.setVel(parcelVelocities);	
				
		if ( verbose_ )
		{							
			int index = 0;
			Info << " Parcel center   " <<  parcelsm.position(index)  << endl;
			Info << " Parcel velocity " <<  parcelsm.velocity(index)  << endl;
			Info << " Parcel diameter " <<  2.*parcelsm.radius(index) << endl;		
		}
		
		// Write parcel diameters in a file --> dparcel_npart=...
		
		// Create output folder	
		if( !isDir(mesh.time().path()/outputRelativePath_) )
		{
			mkDir(mesh.time().path()/outputRelativePath_);
		} 
		
		char charfilterwidth[100];
		sprintf(charfilterwidth,"dparcel_npart=%d",nparticle_);
		fileName outputfile(charfilterwidth);
		OFstream str_dparcel(mesh.time().path()/outputRelativePath_/outputfile);
		str_dparcel << "# Number of parcels = " << parcelI << nl;
		str_dparcel << "# Number of particles not used to construct parcels = " << sm.numberOfParticles()-parcelI*nparticle_ << nl;
		
		//label cellI = 0;
		for(int index = 0; index <  parcelsm.numberOfParticles(); index++)
	    	{
			str_dparcel << 2.*parcelsm.radius(index) << nl;
		
			// Check if filter width is larger? 
			//cellI = parcelsm.cellIDs()[index][0];
			
			/*
			// Calculate filter mesh volume
			const labelList fcellI = stencillistfilter[filterwidth][cellI]; 
			scalar total_volume = 0;											
			forAll(fcellI,filtercellI)
			{
				total_volume += mesh.V()[fcellI[filtercellI]];
			}
			total_volume = pow(total_volume,1./3.);
			
			if ( 2.*parcelsm.radius(index) > total_volume )
			{
				Info << " " 			<< endl;
				Info << " !!! WARNING !!! "	<< endl;
				Info << " Parcel= " << index << " > " << " mesh volume " << total_volume << endl;
				Info << " " 			<< endl;
				//str_dparcel << " " << total_volume << nl;
			}
			*/	
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
		const labelListList& parcelList
	)
 	{
 		// Distance vector
		vector dist(0,0,0);
		
		// Parcel averaged force 
		vectorField parcelForce(parcelsm.numberOfParticles());
		
		// Parcel averaged drag force
		vectorField parcelDragForce(parcelsm.numberOfParticles());
		
		// Parcel averaged filtered relative velocity 
		vectorField parcelFilteredRelativeVelocity(parcelsm.numberOfParticles());
		
		// Parcel averaged relative velocity
		vectorField parcelRelativeVelocity(parcelsm.numberOfParticles());		
		
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
			
		forAll(parcelList,parcelII)	
		{ 			
			labelList parcelListL = parcelList[parcelII];
			dist = vector(0,0,0);

			for(int j=0;j<3;j++) 
			{
						    parcelForce[parcelII][j] = 0; 										
				 parcelFilteredRelativeVelocity[parcelII][j] = 0; 
						parcelDragForce[parcelII][j] = 0;
					 parcelRelativeVelocity[parcelII][j] = 0;
			}					

			forAll(parcelListL,particleI)
			{				
				dist[0] = sqrt(   ( sm.position(parcelListL[particleI]).x() - parcelsm.position(parcelII).x() ) 
					            * ( sm.position(parcelListL[particleI]).x() - parcelsm.position(parcelII).x() ) );
				
				dist[1] = sqrt(   ( sm.position(parcelListL[particleI]).y() - parcelsm.position(parcelII).y() ) 
					            * ( sm.position(parcelListL[particleI]).y() - parcelsm.position(parcelII).y() ) );
				
				dist[2] = sqrt(   ( sm.position(parcelListL[particleI]).z() - parcelsm.position(parcelII).z() ) 
					            * ( sm.position(parcelListL[particleI]).z() - parcelsm.position(parcelII).z() ) );				
		
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
							            parcelForce[parcelII][j] += ( - Volp * gradPg_int[j] 
							                                          + Volp * gradbarPg_int[j]
								 		 	          + DragForce_[parcelListL[particleI]][j]   ) / ( dist[j] + SMALL ) ; 
										
						 parcelFilteredRelativeVelocity[parcelII][j] += Ur[j] / ( dist[j] + SMALL); 
					
						                parcelDragForce[parcelII][j] += DragForce_[parcelListL[particleI]][j] / ( dist[j] + SMALL );
							
						         parcelRelativeVelocity[parcelII][j] += (Ufluid[j] - Up[j]) / ( dist[j] + SMALL ); 										

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
			 parcelbarBetai_[parcelII][j] =        parcelForce[parcelII][j] / ( parcelFilteredRelativeVelocity[parcelII][j] ) ;
			      parcelgii_[parcelII][j] =    parcelbarBetai_[parcelII][j] *           parcelRelativeVelocity[parcelII][j] 
				                       / ( parcelDragForce[parcelII][j] ) ;				
			}
		}
		
		if( verbose_ )
		{
		//Info << parcelForce     << endl;
		//Info << parcelFilteredRelativeVelocity << endl;
		//Info << parcelDragForce << endl;
		//Info << parcelRelativeVelocity << endl;
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
		const scalar& maxalpp_,
		const volScalarField& baralpf_,
		const volVectorField tildeUf_,
		vectorField& gii_,
		Field< Field < Field <scalar> > >& giiCondalppFilter_,
		Field< Field < Field <scalar> > >& numbergiiCondalppFilter_				
	)
	{			
		const scalar minalpp_ = 0.;
		if ( maxalpp_ == 0)
		{
			FatalError << " Max. solid volumde fraction equal to = " << maxalpp_ << abort(FatalError);
		}	
		const scalar binDx_ = (maxalpp_-minalpp_)/nBin_;	
		
		interpolationCellPoint<scalar> baralpfInterpolator_(baralpf_);		
		
		label cellI = 0;
		vector position(0,0,0);
		scalar baralpp = 0;			
		
		for(int index = 0; index <  sm.numberOfParticles(); index++)
	    	{
			cellI = sm.cellIDs()[index][0];
			
			// Info << " cellI " << cellI << endl;
			    
			if (cellI > -1) // particle Found
			{
								
				position = sm.position(index);
				baralpp = 1.-baralpfInterpolator_.interpolate(position,cellI) ;
	    	 		label binI = floor(baralpp/binDx_);			
				
				for(int j=0;j<3;j++) 
				{	
					      giiCondalppFilter_[binI][j][filterwidth] += gii_[index][j];					      
					numbergiiCondalppFilter_[binI][j][filterwidth] += 1;																		
				}
				
			}
		}
		
	} 

	void writeBins
	(
	    	const fvMesh& mesh,
		const int& nBin_,
		const scalar& maxalpp_,
		const fileName outputRelativePath_,		
		const int& minfilterwidth,
		const int& maxfilterwidth,
		const int& filterincrement,
		const int& nparticle_,
		Field< Field < Field <scalar> > >& giiCondalppFilter_,
		Field< Field < Field <scalar> > >& numbergiiCondalppFilter_,
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
		char charfilterwidthN[100]; 
		
		// Debugging for this subroutine
		// verbose_ = true;
						
		for( int filterwidth = minfilterwidth; filterwidth <= maxfilterwidth; filterwidth += filterincrement )
		{

			int Filter = 2*(filterwidth+1)+1;

			sprintf(charfilterwidth, "h_vs_baralpp_%dX%dX%d",Filter,Filter,Filter);
			fileName outputfile(charfilterwidth);
				
			sprintf(charfilterwidthN, "Nrealization_%dX%dX%d",Filter,Filter,Filter);
			fileName outputfileN(charfilterwidthN);
							
			if ( nparticle_ > 0 )
			{								
				sprintf(charfilterwidth,"h_vs_baralpp_%dX%dX%d_npart=%d",Filter,Filter,Filter,nparticle_);
				outputfile = charfilterwidth;
				sprintf(charfilterwidthN,"Nrealization_%dX%dX%d_npart=%d",Filter,Filter,Filter,nparticle_);
				outputfileN = charfilterwidthN;				
			}
			
			OFstream str_gii(mesh.time().path()/outputRelativePath_/outputfile);
			OFstream str_nreal(mesh.time().path()/outputRelativePath_/outputfileN);
			
			Info << " 		   " << endl; 		
			Info << " Writing the file " << outputfile << ", " << outputfileN 
			     << " into the folder: " << mesh.time().path()/outputRelativePath_ << endl;
			
			str_gii   << "# baralpp \t H_x \t H_y \t H_z"  << nl;
			str_nreal << "# baralpp \t Nx  \t Ny  \t H_z"  << nl;
			
			for( int i = 0; i < nBin_; i++ )
			{


				str_nreal << (i+1./2.)*maxalpp_/nBin_  << " " << numbergiiCondalppFilter_[i][0][filterwidth] 
				                        	       << " " << numbergiiCondalppFilter_[i][1][filterwidth] 
								       << " " << numbergiiCondalppFilter_[i][2][filterwidth]		<< nl;
								       
								       
				if ( verbose_ )				       
				{
					Info  <<  (i+1./2.)*maxalpp_/nBin_  << " " << giiCondalppFilter_[i][0][filterwidth] << " " << numbergiiCondalppFilter_[i][0][filterwidth] 
							        	    << " " << giiCondalppFilter_[i][1][filterwidth] << " " << numbergiiCondalppFilter_[i][1][filterwidth] 
							        	    << " " << giiCondalppFilter_[i][2][filterwidth] << " " << numbergiiCondalppFilter_[i][2][filterwidth] 	<< endl;
				}
				
				// Avoid to divide by zero
				for(int j=0;j<3;j++) 
				{				  						
					if ( numbergiiCondalppFilter_[i][j][filterwidth] == 0 ) 
					{
						numbergiiCondalppFilter_[i][j][filterwidth] = 1; 
 						      giiCondalppFilter_[i][j][filterwidth] = 0;
					}
				}					       		

				str_gii  <<  (i+1./2.)*maxalpp_/nBin_  << " " << giiCondalppFilter_[i][0][filterwidth] / numbergiiCondalppFilter_[i][0][filterwidth] 
					 	   		       << " " << giiCondalppFilter_[i][1][filterwidth] / numbergiiCondalppFilter_[i][1][filterwidth] 
								       << " " << giiCondalppFilter_[i][2][filterwidth] / numbergiiCondalppFilter_[i][2][filterwidth] 	<< nl;	

				if ( verbose_ )				       
				{
					Info  <<  (i+1./2.)*maxalpp_/nBin_  << " " << giiCondalppFilter_[i][0][filterwidth] / numbergiiCondalppFilter_[i][0][filterwidth] 
							        	    << " " << giiCondalppFilter_[i][1][filterwidth] / numbergiiCondalppFilter_[i][1][filterwidth] 
							        	    << " " << giiCondalppFilter_[i][2][filterwidth] / numbergiiCondalppFilter_[i][2][filterwidth] 	<< endl;

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
	
	
	// Create parcel cloud	
	cfdemCloud parcelCloud(mesh);
	parcelCloud.reAllocArrays();	
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
    	fileName outputRelativePath(propsDict.lookup("outputRelativePath"));
	fileName stencilsOutputRelativePath(propsDict.lookup("stencilsOutputRelativePath"));
	bool particleCoarsening(false);
	bool verboseParticleCoarsening(false);
	//int npartParticleCoarsening(1);
	labelList npartParticleCoarsening(1);
	if(propsDict.found("verbose")) verbose = true;	
	if(propsDict.found("ParticleCoarsening")) 
	{
		particleCoarsening = true;
		if(propsDict.found("verboseParticleCoarsening")) verboseParticleCoarsening = true; 
		npartParticleCoarsening = propsDict.lookup("npartParticleCoarsening");
		labelList test(propsDict.lookup("npartParticleCoarsening"));
	}
	
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
			Foam::DragForce
			(
				particleCloud,
				voidfraction,
				U,
				rho,
				verbose,
				DragForce
			);

			if(particleCoarsening)
			{
				forAll(npartParticleCoarsening,npartI)
				{
					Info << " " << endl;
					Info << " Particle coarsening: "<< endl;

					// Parcel particle list
					labelList particleID(npartParticleCoarsening[npartI],(0));
					labelListList parcelParticleList(particleCloud.numberOfParticles(),particleID);
					
					Foam::createParcels
					(
						mesh,
						particleCloud,
						parcelCloud,
						npartParticleCoarsening[npartI],
						verboseParticleCoarsening,
						outputRelativePath,
						parcelParticleList					
					);
					

					for(int FilterWidth = minFilterWidth; FilterWidth <= maxFilterWidth; FilterWidth += FilterIncrement )
					{
						Info << " " << endl;
						int Filter = 2*(FilterWidth+1)+1;
						Info << " Filter size = " << Filter << "X" << Filter << "X" << Filter << endl;

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
							barPg
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

						// Particle coarsening active or not
						Foam::ParticleCoarsening
						(
							mesh,
							StencilListFilter,
							FilterWidth,
							particleCloud,
							parcelCloud,
							verboseParticleCoarsening,
							U,
							tildeUf,
							p,
							barPg,
							DragForce,
							ResolvedDragForce,
							parcelgii,
							outputRelativePath,
							parcelParticleList					
						);

						Foam::applyBins
						(				
							parcelCloud,
							FilterWidth,
							maxFilterWidth,
							nBin,
							maxalpp,
							baralpf,
							tildeUf,
							parcelgii,
							NparcelgiiCondalppFilter[npartI],
							NparcelnumbergiiCondalppFilter[npartI]	
						);				

						Foam::applyBins
						(				
							particleCloud,
							FilterWidth,
							maxFilterWidth,
							nBin,
							maxalpp,
							baralpf,
							tildeUf,
							gii,
							giiCondalppFilter,
							numbergiiCondalppFilter
						);					

					}
				}	
			
			}
			else
			{
				for(int FilterWidth = minFilterWidth; FilterWidth <= maxFilterWidth; FilterWidth += FilterIncrement )
				{
					Info << " " << endl;
					int Filter = 2*(FilterWidth+1)+1;
					Info << " Filter size = " << Filter << "X" << Filter << "X" << Filter << endl;

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
						barPg
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

					Foam::applyBins
					(				
						particleCloud,
						FilterWidth,
						maxFilterWidth,
						nBin,
						maxalpp,
						baralpf,
						tildeUf,
						gii,
						giiCondalppFilter,
						numbergiiCondalppFilter
					);
				}		
			}	
		
		}	

	}
	
	// Write bins		
	Foam::writeBins
	(
		mesh,
		nBin,
		maxalpp,		
		outputRelativePath,
		minFilterWidth, 
		maxFilterWidth,
		FilterIncrement,
		-1,				// Just for fluid coarsening results
		giiCondalppFilter,
		numbergiiCondalppFilter,
		verbose						
	);
	
	if(particleCoarsening)
	{
		forAll(npartParticleCoarsening,npartI)
		{
			Foam::writeBins
			(
				mesh,
				nBin,
				maxalpp,		
				outputRelativePath,
				minFilterWidth, 
				maxFilterWidth,
				FilterIncrement,
				npartParticleCoarsening[npartI],
				NparcelgiiCondalppFilter[npartI],
				NparcelnumbergiiCondalppFilter[npartI],
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
