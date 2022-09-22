/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Copyright of this Contribution
    Copyright 2014-     Ali Ozel & Sankaran Sundaresan, Princeton University

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "OzelSundaresanDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(OzelSundaresanDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    OzelSundaresanDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
OzelSundaresanDrag::OzelSundaresanDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    granVelFieldName_(propsDict_.lookup("granVelFieldName")),
    Us_(sm.mesh().lookupObject<volVectorField> (granVelFieldName_)),    
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
//_AO 06/02/2014 
//	phi_(readScalar(propsDict_.lookup("phi"))),
//_AO 06/02/2014 
    interpolation_(false),
    scaleDia_(1.),
    scaleDrag_(1.),
    switchingVoidfraction_(0.8),
    velslipdiff_(false),
    ErgunMix_(false),     
    Kxx_(3.33),
    Kyy_(3.33),
    Kzz_(3.33),    
    const_h1(0.2),
    const_h2(-1.88),
    const_h3(5.16),
    alpsmax(0.64),
        gravityFieldName_(propsDict_.lookup("gravityFieldName")),
    #if defined(version21) || defined(version16ext)
        g_(sm.mesh().lookupObject<uniformDimensionedVectorField> (gravityFieldName_))
    #elif defined(version15)
        g_(dimensionedVector(sm.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(gravityFieldName_)).value())
    #endif
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "OzelSundaresanDrag.logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must  be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().scalarFields_.append("Rep");
    particleCloud_.probeM().scalarFields_.append("betaP");
    particleCloud_.probeM().scalarFields_.append("voidfraction");
    particleCloud_.probeM().writeHeader();

//_AO 06/02/2014
    //if (propsDict_.found("verbose")) verbose_=readBool(propsDict_.lookup("verbose"));
    if (propsDict_.found("verbose")) verbose_=true;
//_AO 06/02/2014       
   
    if (propsDict_.found("treatExplicit")) treatExplicit_=readBool(propsDict_.lookup("treatExplicit"));
    if (propsDict_.found("interpolation")) interpolation_=readBool(propsDict_.lookup("interpolation"));
    if (propsDict_.found("implDEM"))
    {
        implDEM_=readBool(propsDict_.lookup("implDEM"));
        treatExplicit_=false;
        setImpDEMdrag();
        Info << "Using implicit DEM drag formulation." << endl;
    }
    particleCloud_.checkCG(true);
    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_=scalar(readScalar(propsDict_.lookup("scaleDrag")));

    if (propsDict_.found("switchingVoidfraction"))
        switchingVoidfraction_ = readScalar(propsDict_.lookup("switchingVoidfraction"));

    if (propsDict_.found("velslipdiff")) 
    {
       velslipdiff_     = readBool(propsDict_.lookup("velslipdiff"));
       interpolation_=true; // Force calculate at the particle position
       Info << "OzelSundaresanDrag: Calculation of local relative velocity at grid location by using Eulerian velocities then interpolate to particle position is ACTIVE." << endl
              << "OzelSundaresanDrag: Also, interpolation at the particle position hase been FORCED TO BE ACTIVE." << endl;
    }  

    if (propsDict_.found("Kxx"))
        Kxx_=readScalar(propsDict_.lookup("Kxx"));
    if (propsDict_.found("Kxx"))
        Kyy_=readScalar(propsDict_.lookup("Kyy"));
    if (propsDict_.found("Kzz"))
        Kzz_=readScalar(propsDict_.lookup("Kzz"));    
    if (propsDict_.found("const_h1"))
        const_h1=readScalar(propsDict_.lookup("const_h1"));


    Info << "OzelSundaresan - interpolation switch: " << interpolation_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

OzelSundaresanDrag::~OzelSundaresanDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void OzelSundaresanDrag::setForce() const
{
    if (scaleDia_ > 1)
    {
        Info << "OzelSundaresan using scale = " << scaleDia_ << endl;
//_AO 06/02/2014 
	//scaleDrag_ = scaleDia_ * scaleDia_ * scaleDia_ ; 
	//Info << "OzelSundaresan using drag scale = " << scaleDrag_ << endl;
//_AO 06/02/2014 
    }
    else if (particleCloud_.cg() > 1){
        scaleDia_=particleCloud_.cg();
        Info << "OzelSundaresan using scale from liggghts cg = " << scaleDia_ << endl;
//_AO 06/02/2014 
	//scaleDrag_ = scaleDia_ * scaleDia_ * scaleDia_ ; 
	//Info << "OzelSundaresan using drag scale = " << scaleDrag_ << endl;
//_AO 06/02/2014 
    }

    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu() / rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI=0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Vs(0);
    scalar localPhiP(0);

    scalar CdMagUrLag(0);       //Cd of the very particle
    scalar KslLag(0);           //momentum exchange of the very particle (per unit volume)

    // Kronecker delta and particle resistance tensor
    tensor deltaij(1,  0,  0,
                   0,  1,  0,
                   0,  0,  1);
    tensor Kii = tensor(Kxx_,    0,    0,     
                           0, Kyy_,    0,
                           0,    0, Kzz_);
    // Drag coefficient tensor
    scalar betaP(0.0);
    tensor betaPTensorial(0,0,0,
                          0,0,0, 
                          0,0,0);
    //temporary drag correction factors
    scalar func_f(1.0);
    scalar func_h(1.0); 

   interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    // Slip velocity 
    interpolationCellPoint<vector> UsInterpolator_(Us_);	

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            
            cellI = particleCloud_.cellIDs()[index][0];
           
            //reset forces and force tensors
            drag = vector(0,0,0);
    	    tensor Kii = tensor(Kxx_,    0,    0, 
	                           0, Kyy_,    0,
	                           0,    0, Kzz_);
            betaP = 0.0;
            betaPTensorial = tensor(0,0,0,
                                    0,0,0, 
                                    0,0,0);
            Vs = 0;
            Ufluid =vector(0,0,0);
            voidfraction=0;

            if (cellI > -1) // particle Found
            {

                if(interpolation_)
                {
	            position     = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid       = UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;
                    if(voidfraction>1.00) voidfraction = 1.0;
                    if(voidfraction<0.10) voidfraction = 0.10;
                }
                else
                {
		    voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                }

                //slit velocity calculation

        	if(velslipdiff_)
                {
    		    	Ur =  Ufluid
                          - UsInterpolator_.interpolate(position,cellI);	
                }
                else
                {
                    Us = particleCloud_.velocity(index);
                    Ur = Ufluid - Us;
                }

                magUr = mag(Ur);
                ds = 2*particleCloud_.radius(index); //_AO 06/02/2014 *phi_; Why ???
                rho = rho_[cellI];
                nuf = nufField[cellI];

                Rep=0.0;
                localPhiP = 1.0f-voidfraction+SMALL;
                Vs = ds*ds*ds*M_PI/6;

                //Compute specific drag coefficient (i.e., Force per unit slip velocity and per m³ SUSPENSION)
                // if(voidfraction > switchingVoidfraction_) //dilute
                //{
                    Rep=ds/scaleDia_*voidfraction*magUr/nuf;
                    CdMagUrLag = (24.0*nuf/(ds/scaleDia_*voidfraction)) //1/magUr missing here, but compensated in expression for KslLag!
                                 *(scalar(1)+0.15*Foam::pow(Rep, 0.687));

                    KslLag = 0.75*(
                                            rho*localPhiP*voidfraction*CdMagUrLag
                                          /
                                            (ds/scaleDia_*Foam::pow(voidfraction,2.65))
                                          );
                //}
		/*
                else if( voidfraction < switchingVoidfraction_ && ErgunMix_ ) //dense
                {
                    KslLag = (150*localPhiP*localPhiP*nuf*rho)/
                             (voidfraction*ds/scaleDia_*ds/scaleDia_+SMALL)
                            +
                             (1.75*(localPhiP) * magUr * rho)/
                             ((ds/scaleDia_));
                }
		*/
    		    //Compute a Dimensionless characteristic mesh length and mesh TODO
//                scalar taup = 0.75 * rhos/rho
//                         *(ds/scaleDia_) 
//                          /(CD*magUr + SMALL);

     	      // Dimensionless characteristic mesh length
//		      delta_bar = pow(particleCloud_.mesh().V()[cellI],1./3.)/(taup*magUr);
//		      func_f = ???; //A posteriori test Delta_f=9X9X9 //pow(delta_bar,2)/(const_f1+pow(delta_bar,2));	TODO 		
		     
                 // Volume fraction dependency of drag correction
 		         func_h = - tanh(localPhiP/const_h1)
                             * sqrt(localPhiP/alpsmax)*(1.0-localPhiP/alpsmax)*(1.0-localPhiP/alpsmax)
                             * (    
                                    1.0
                                   +const_h2    * localPhiP/alpsmax
                                   +const_h3    *(localPhiP/alpsmax)*(localPhiP/alpsmax)
                                );

		// Limiting Kii*func_h*func_f <=-0.999
		if( Kii.xx()*func_h*func_f <=-0.999 )
    		       Kii.xx() = -0.999/(func_h*func_f);
 
                if( Kii.yy()*func_h*func_f <=-0.999 )
                       Kii.yy() = -0.999/(func_h*func_f);

                if( Kii.zz()*func_h*func_f <=-0.999 )
                       Kii.zz() = -0.999/(func_h*func_f);

                // calc particle's drag coefficient (i.e., Force per unit slip velocity and per m³ PARTICLE)
                betaP = KslLag / localPhiP;
                betaPTensorial = betaP  
                                        * (
                                              deltaij                    //isotropic part
                                           + Kii*func_h*func_f  //off-diagonal part
                                           ) ;

                // calc particle's drag, Vs is the total volume of the parcel
                drag = Vs * (betaPTensorial & Ur) * scaleDrag_;

                if (modelType_=="B")
                    drag /= voidfraction;

                if(verbose_ && index >=0 && index <2)
                {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "ds/scale = " << ds/scaleDia_ << endl;
                    //_AO 06/02/2014 
		    // Pout << "phi = " << phi_ << endl;
                    //_AO 06/02/2014 
		    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
        	    Pout << "func_f = " << func_f << endl;
                    Pout << "func_h = " << func_h << endl;
                    Pout << "Kii = " << Kii << endl;		
                    Pout << "betaP = " << betaP << endl;
                    Pout << "betaPTensorial = " << betaPTensorial << endl;
                    Pout << "drag = " << drag << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(drag);   //first entry must the be the force
                    vValues.append(Ur);
                    sValues.append(Rep);
                    sValues.append(betaP);
                    sValues.append(voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // set force on particle
            if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += drag[j];
            else  for(int j=0;j<3;j++) impForces()[index][j] += drag[j];

            // set Cd
            if(implDEM_)
            {
                for(int j=0;j<3;j++) fluidVel()[index][j]=Ufluid[j];

                if (modelType_=="B" && cellI > -1)
                    Cds()[index][0] = Vs*betaP/voidfraction*scaleDrag_;
                else
                    Cds()[index][0] = Vs*betaP*scaleDrag_;

            }else{
                for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];
            }

        //}// end if mask
    }// end loop particles
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
