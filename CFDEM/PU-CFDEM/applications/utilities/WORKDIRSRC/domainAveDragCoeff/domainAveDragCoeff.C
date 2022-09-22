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

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"

#include "cfdemCloud.H"
#include "dataExchangeModel.H"
#include "voidFractionModel.H"
#include "locateModel.H"
#include "averagingModel.H"
#include "momCoupleModel.H"
#include "forceModel.H"
#include "IOModel.H"
#include "interpolationCellPoint.H"

#include "timeSelector.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

        	if(ds_ <= 0) Pout << " ds " << ds_ << endl;
		if(voidfraction_> 1.0 ) Pout << " Int. voidfraction > 1 " << " value= " << voidfraction_ << endl;				
		
		if (magUr > 0 && ds_ > 0)
		{
			// calc particle Re number
			Rep = voidfraction_*ds_*magUr/(nuf_+SMALL);
            
			// calc CD
			if (Rep < 1000)
			{
				CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction_,-1.65);
			}
			else
			{
				CD = 0.44*pow(voidfraction_,-1.65); 
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

		
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{	
    	timeSelector::addOptions();
	#   include "addRegionOption.H"
	
    	argList::validArgs.append("rhop");
    	argList::validArgs.append("driftVel");
		
	#   include "setRootCase.H"
	#   include "createTime.H"
    	instantList timeDirs = timeSelector::select0(runTime, args);
	#   include "createNamedMesh.H"
	
	//argList args(argc,argv);
	const scalar rhop = args.argRead<scalar>(1);
   	const bool driftVel = args.argRead<bool>(2);
	
   	// read gravity
    	#include "readGravitationalAcceleration.H"
	
	// Create field variables
	#include "createFields.H" 	 

	// Create particle cloud			
	cfdemCloud particleCloud(mesh);
	
	particleCloud.reAllocArrays();

	double **positions;
	double **velocities;
	double **radii;
	double **cellID;
	
	double **forces;

	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
	
	particleCloud.dataExchangeM().allocateArray(forces,0.,3);	

	#ifdef comp
	    const volScalarField nufField = particleCloud.turbulence().mu()/rho_;
	#else
	    const volScalarField& nufField = particleCloud.turbulence().nu();
	#endif
	
	// Debugging
	bool verbose(true);	

	// Open file
	OFstream* sPtrBeta;
	sPtrBeta =  new OFstream("filteredBeta");
    	*sPtrBeta  << "#Time \t forceBalanceBeta \t domainAveBeta \t vSlipEuler" << endl;		

	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();
		
		// Read gas velocity 
		IOobject Uheader
		(
		   "U",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		   
		Info<< " 	Reading U" << endl;
		volVectorField U(Uheader,mesh);

		// Read gas volume fraction
		IOobject voidfractionHeader
		(
		   "voidfraction",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		   
		Info<< " 	Reading gas volume fraction" << endl;
		volScalarField alpf(voidfractionHeader,mesh);		

		// Read gas pressure gradient
		IOobject pHeader
		(
		   "p",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		   
		Info<< " 	Reading Pg" << endl;
		volScalarField p(pHeader,mesh);

		// Read gas pressure gradient
		IOobject rhoHeader
		(
		   "rho",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);
		   
		Info<< " 	Reading gas density" << endl;
		volScalarField rho(rhoHeader,mesh);
		
		// Domain average gas velocity
		vector domainAveUf(0,0,0);
		domainAveUf = fvc::domainIntegrate(alpf*U).value()/fvc::domainIntegrate(alpf).value();

		// Domain average gas pressure
		volVectorField gradPg(fvc::grad(p));
		vector domainAveGradPg(0,0,0);
		domainAveGradPg = fvc::domainIntegrate(alpf*gradPg).value()/fvc::domainIntegrate(alpf).value();
		interpolationCellPoint<vector> gradPgInterpolator(gradPg);
		
		scalar domainAvePg(0);
		domainAvePg = fvc::domainIntegrate(alpf*p).value()/fvc::domainIntegrate(alpf).value();
		
		// Domain average solid volume fraction
		scalar domainAveAlpf(0);
		domainAveAlpf = fvc::domainIntegrate(alpf).value()/sum(mesh.V()).value();
		
		Info << "\n Mean gas velocity      = " << domainAveUf   << endl;
		Info << " Mean pressure gradient   = " << domainAveGradPg   << endl;
		Info << " Mean gas pressure        = " << domainAvePg   << endl;		
		Info << " Mean gas volume fraction = " << domainAveAlpf << endl;
		
		// Drift velocity
		vector domainAveDrift(0,0,0);
		domainAveDrift = fvc::domainIntegrate((1.-alpf)*U).value()/fvc::domainIntegrate(1.-alpf).value() - domainAveUf;
		
		Info << " Drift velocity = " << domainAveDrift << endl;
		
		if ( runTime.timeName() != "0" )
		{
			
			int count = runTime.value() / particleCloud.dataExchangeM().DEMts();				

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

			particleCloud.dataExchangeM().getData("f","vector-atom",forces,count);
			Info<< " 	Reading forces on particles" << endl;
			Info<< " " << endl;				
		
			// Locate particles in Eulerian grid
			particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
			particleCloud.setPos(positions);
			particleCloud.setVel(velocities);

			if(verbose)
			{
				labelList exIndex(1,1);
				forAll(exIndex,ii)
				{
					int index = exIndex[ii];
					
					Info << tab << "" << endl;
					Info << tab << "index  = " << index << endl;
					Info << tab << "rp     = " << particleCloud.radius(index) << endl;
					Info << tab << "Vp     = " << particleCloud.velocity(index) << endl;
					Info << tab << "Xp     = " << particleCloud.position(index) << endl;
					Info << tab << "CellID = " << particleCloud.particleCell(index) << endl;
					Info << tab << "Fp     = " << forces[index][0] << " " << forces[index][1] << " " << forces[index][2] << endl;
					if ( particleCloud.particleCell(index) > -1 )
						Info << " Uf     = " << U[particleCloud.particleCell(index)] << endl;
					Info << "" << endl;
				}	
			}
			
			label cellI(-1);
			vector position(0,0,0);
			vector gradPgXp(0,0,0);
			vector gradPg(0,0,0);
			vector Ufluid(0,0,0);
			vector Up(0,0,0);
			vector Ur(0,0,0);
			scalar ds(0);
			scalar Volp(0);
			scalar voidfraction(0);
			scalar nuf(0);
			scalar rhof(0);
			vector totalUp(0,0,0);
			vector WenYuDrag(0,0,0);
			vectorField parFilteredBeta(particleCloud.numberOfParticles(),vector(0,0,0));
			vector totalParFilteredBeta(0,0,0);
			label nP(0);
			label nPIn(0);
			scalar nPOutOfDomain(0);
			
			for(int parI=0; parI < particleCloud.numberOfParticles(); parI++ )
			{
           
        	    		// Main particle/parcel variables
        	    		cellI = particleCloud.cellIDs()[parI][0];
		    			
        	    		if( cellI > -1 )
               	    		{			    
			    		nPIn++;
			    		  position = particleCloud.position(parI);
					  gradPgXp = gradPgInterpolator.interpolate(position,cellI); 
					    gradPg = domainAveGradPg;
        		    		    Ufluid = U[cellI]; 
        		    		        Up = particleCloud.velocity(parI);
			    			Ur = Ufluid - Up;
        		    			ds = 2*particleCloud.radius(parI);			    			    
        		                      Volp = ds*ds*ds*M_PI/6;
				      voidfraction = alpf[cellI]; 
					       nuf = nufField[cellI];
        		                      rhof = rho[cellI];
					      
					   totalUp+= particleCloud.velocity(parI);
					
					// Calculate drag force
					WenYuDragForce(Ur,ds,rhof,nuf,voidfraction,WenYuDrag);			    

					for(int j=0;j<3;j++)
                			{
                		    	   
					   parFilteredBeta[parI][j] =   ( Volp * gradPgXp[j] - Volp * domainAveGradPg[j] + WenYuDrag[j] ) 
					                              / ( domainAveUf[j] - Up[j] + SMALL );
	
					   if( (parFilteredBeta[parI][j] < scalar(1+SMALL) ) && (parFilteredBeta[parI][j] > scalar(-SMALL) ) ) 
					    {
					    	totalParFilteredBeta[j]+= parFilteredBeta[parI][j];
					    	nP++;
					    }	
					}
				}else
				{
					nPOutOfDomain++;
					//Info << " Particle " << parI << " is not in the domain " << endl;
				}				
			
			}
			
			Info << " Number of particles out of domain = " << nPOutOfDomain << endl; 
			
			vector filteredBetaDomAveVar(0,0,0);
			filteredBetaDomAveVar = totalParFilteredBeta / nP; 
			Info << "\nAveraged filtered drag coefficient = " << filteredBetaDomAveVar << "\n" <<endl;
			
			Info << "\nCalculating filtered drag coefficient by force balance ... " << endl;
			vector vSlip(0,0,0);
			vSlip = domainAveUf - totalUp/nPIn; 
			if(driftVel) vSlip += domainAveDrift; 
			Info << " Domain-ave slip = " << vSlip << endl;
			Info << " Mean particle velocity = " << totalUp/nPIn << endl;
    			Info << " Particle density = " << rhop << endl;
			
			vector filteredBetaForceCalc(0,0,0);
			for(int j=0;j<3;j++) filteredBetaForceCalc[j] = domainAveAlpf * (rhop-rhof) * mag(g.value()) * M_PI/6. * ds*ds*ds / vSlip[j];
			Info << " Filtered drag coefficient = " << filteredBetaForceCalc << endl;
			
			Info << "\nCalculating coarse beta ... " << endl;
			vector domainAveVr(0,0,0);
			domainAveVr = domainAveUf - totalUp/nPIn; 
			vector WenYuDragCoarse(0,0,0);
			WenYuDragForce(domainAveVr,ds,rhof,nuf,domainAveAlpf,WenYuDragCoarse);
			vector coarseBeta(0,0,0);
			for(int j=0;j<3;j++) coarseBeta[j] = WenYuDragCoarse[j] / domainAveVr[j]; 
			Info << " Coarse beta = " << coarseBeta << endl; 
			
			Info << "\nNormalized filtered drag coefficient by " << endl;
			Info << " Force balance calculation " << filteredBetaForceCalc[0] / coarseBeta[0] << " " << filteredBetaForceCalc[1] / coarseBeta[1] << " " << filteredBetaForceCalc[2] / coarseBeta[2] << endl;
			Info << " Domain-averaged values    " << filteredBetaDomAveVar[0] / coarseBeta[0] << " " << filteredBetaDomAveVar[1] / coarseBeta[1] << " " << filteredBetaDomAveVar[2] / coarseBeta[2] << endl;
			
			// Write into file
			    *sPtrBeta << setw(IOstream::defaultPrecision() + 6) 
                		      << U.mesh().time().value()     			<< "   " 
                		      << filteredBetaForceCalc[0] / coarseBeta[0] 	<< "   "
                		      << filteredBetaForceCalc[1] / coarseBeta[1] 	<< "   "
                		      << filteredBetaForceCalc[2] / coarseBeta[2] 	<< "   "		
                		      << filteredBetaDomAveVar[0] / coarseBeta[0] 	<< "   "
                		      << filteredBetaDomAveVar[1] / coarseBeta[1] 	<< "   "
                		      << filteredBetaDomAveVar[2] / coarseBeta[2] 	<< "   "
				      << vSlip[0] << " " << vSlip[1] << " " << vSlip[2] << endl;
		}	


	}
		
	particleCloud.dataExchangeM().destroy(positions,3);
	particleCloud.dataExchangeM().destroy(velocities,3);

	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
	
	

}


// ************************************************************************* //
