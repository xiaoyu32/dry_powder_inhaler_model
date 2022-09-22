#include "domainFilter.H"

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

		Pout << "Domain " << " x=[" << min_x << ":" << max_x << "]" << " y=[" << min_y << ":" << max_y << "]" << " z=[" << min_z << ":" << max_z << "]" << endl;

		const scalar cell_volume_ref = mesh.V()[0];
		forAll(mesh.cells(),cellI)
		{
			if ( abs(cell_volume_ref-mesh.V()[cellI]) > SMALL)  FatalError<< "Non-uniform mesh !!! " << abort(FatalError);
		}

		const scalar delta = pow(cell_volume_ref,1./3);
		Pout << "Uniform mesh,  volume of a cell = " << cell_volume_ref << endl;
		Pout << "Delta_x = Delta_y = Delta_z = " << delta << endl;


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
		Pout << "Number of coupled patches = " << nCoupledPatches << endl;
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

}	