#include "readVariables.H"

namespace Foam
{	
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
	
	void readFilteredEulerianVariables
	(
	    const argList& args, 
	    const Time& runTime, 
	    const fvMesh& mesh,
		volScalarField& alpf_,
		volVectorField& Uf_,
		volScalarField& p_,
		const int& filterwidth_	
	)
	{		
		char charfilterwidth[100]; 
		int Filter = 2*(filterwidth_+1)+1;
		
		// Read filtered gas volume fraction		
		sprintf(charfilterwidth, "barvoidfraction_%dX%dX%d",Filter,Filter,Filter);
		IOobject voidfractionheader
		(
		   charfilterwidth,
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		Pout<< " 	Reading " << charfilterwidth << endl;
		volScalarField voidfraction_(voidfractionheader,mesh);	
		alpf_ = voidfraction_ ;		
		
		
		// Read filtered gas velocity 
		sprintf(charfilterwidth, "tildeU_%dX%dX%d",Filter,Filter,Filter);	
		IOobject Uheader
		(
		   charfilterwidth,
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);	   
		Pout<< " 	Reading " << charfilterwidth << endl;
		volVectorField U_(Uheader,mesh);
		Uf_ = U_ ;
			
		
		// Read filtered gas pressure
		sprintf(charfilterwidth, "barPg_%dX%dX%d",Filter,Filter,Filter);
		IOobject pheader
		(
		   charfilterwidth,
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		Pout<< " 	Reading " << charfilterwidth << endl;
		volScalarField Pg_(pheader,mesh);
		p_ = Pg_ ;	
		
	}	
	
}	
