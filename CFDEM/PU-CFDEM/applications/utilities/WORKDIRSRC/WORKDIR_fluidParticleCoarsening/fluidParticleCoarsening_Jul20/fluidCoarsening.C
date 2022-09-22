#include "fluidCoarsening.H"

namespace Foam
{

	void filteringEulerianScalar
	(
		const fvMesh& mesh,
		const labelListListList& stencillistfilter,
		const int& filterwidth,
		const volScalarField& voidfraction_,
		const volScalarField& phi_,
		volScalarField& fphi_ 
	)
	{

		forAll(mesh.cells(),cellI)
		{
			scalar barvoidfraction = 0;
			fphi_[cellI] = 0; 			
			const labelList fcellI = stencillistfilter[filterwidth][cellI]; 

			forAll(fcellI,filtercellI)
			{
				barvoidfraction +=  voidfraction_[fcellI[filtercellI]] 
					               * mesh.V()[fcellI[filtercellI]]; 				
				fphi_[cellI]    +=  voidfraction_[fcellI[filtercellI]] 
					               *     phi_[fcellI[filtercellI]] 
						       * mesh.V()[fcellI[filtercellI]]; 
			}
			if ( barvoidfraction > 0 )
			{
				fphi_[cellI] /= barvoidfraction;
			}
			else 
			{
				fphi_[cellI] = 0;
			}
		} 

	}	

	void filteringEulerianVector
	(
		const fvMesh& mesh,
		const labelListListList& stencillistfilter,
		const int& filterwidth,
		const volScalarField& voidfraction_,
		const volVectorField& phi_,
		volVectorField& fphi_ 
	)
	{

		forAll(mesh.cells(),cellI)
		{
			scalar barvoidfraction = 0;
			fphi_[cellI] = vector(0,0,0); 			
			const labelList fcellI = stencillistfilter[filterwidth][cellI]; 

			forAll(fcellI,filtercellI)
			{
				barvoidfraction +=  voidfraction_[fcellI[filtercellI]] 
					               * mesh.V()[fcellI[filtercellI]]; 				
				fphi_[cellI]    +=  voidfraction_[fcellI[filtercellI]] 
					               *     phi_[fcellI[filtercellI]] 
						       * mesh.V()[fcellI[filtercellI]]; 
			}
			if ( barvoidfraction > 0 )
			{
				fphi_[cellI] /= barvoidfraction;
			}
			else 
			{
				fphi_[cellI] = vector(0,0,0);
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
		const bool EulerEulerFiltering_,
		const bool EulerianVslipBin_  
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
		
		filteringEulerianVector
		(
			mesh,
			stencillistfilter,
			filterwidth,
			voidfraction_,
			U_,
			tildeU 
		);
			
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
		
		filteringEulerianScalar
		(
			mesh,
			stencillistfilter,
			filterwidth,
			voidfraction_,
			p_,
			barp 
		);
				
		Info<< " 	Writing filtered gas pressure" << endl;
	    	barp.write();		
		barPg_ = barp;	
				
		if( EulerEulerFiltering_ || EulerianVslipBin_ )
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
			
			filteringEulerianVector
			(
				mesh,
				stencillistfilter,
				filterwidth,
				1-voidfraction_,
				Up_,
				tildeUp 
			);

			Info<< " 	Writing filtered Up" << endl;
	    		tildeUp.write();		
			tildeUp_ = tildeUp;		
		}

    }
}	