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
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "WenYuANNDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"
#include <fstream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WenYuANNDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    WenYuANNDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
WenYuANNDrag::WenYuANNDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    multiTypes_(false),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
//    phi_(propsDict_.lookupOrDefault<scalar>("phi",1.)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    maxTypeCG_(1),
    typeCG_(propsDict_.lookupOrDefault<scalarList>("coarseGrainingFactors",scalarList(1,1.0))),
    scaleDia_(1.),
    scaleDrag_(1.)
    //,switchingVoidfraction_(propsDict_.lookupOrDefault("switchingVoidfraction", 0.9))
{
    //Append the field names to be probed
    // suppress particle probe
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");
        particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must  be the force
        particleCloud_.probeM().vectorFields_.append("Urel");
        particleCloud_.probeM().scalarFields_.append("Rep");
        particleCloud_.probeM().scalarFields_.append("betaP");
        particleCloud_.probeM().scalarFields_.append("voidfraction");
        particleCloud_.probeM().writeHeader();
    }

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
/*    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_IMPL_FORCE_DEM,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(SW_SCALAR_VISCOSITY,true); // activate scalarViscosity switch
*/



    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(2,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch


 
    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(true);
    if (propsDict_.found("scale") && typeCG_.size()==1)
    {
        // if "scale" is specified and there's only one single type, use "scale"
        scaleDia_ = scalar(readScalar(propsDict_.lookup("scale")));
        typeCG_[0] = scaleDia_;
    }
    if (propsDict_.found("scaleDrag"))
    {
        scaleDrag_=scalar(readScalar(propsDict_.lookup("scaleDrag")));
    }

    if (typeCG_.size() > 1)
    {
        multiTypes_ = true;
    }

    //-------------------Added for NN--------------//
    if (propsDict_.found("useReNNCorr")) {
        nnDragCorrectorPtr_.reset(new NNDragCorrector(propsDict_, sm.mesh()));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

WenYuANNDrag::~WenYuANNDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WenYuANNDrag::setForce() const
{
    if (typeCG_.size() > 1 || typeCG_[0] > 1)
    {
        Info<< "WenYuANN using scale = " << typeCG_ << endl;
    }
    else if (particleCloud_.cg() > 1.)
    {
        scaleDia_ = particleCloud_.cg();
        typeCG_[0] = scaleDia_;
        Info<< "WenYuANN using scale from liggghts cg = " << scaleDia_ << endl;
    }

    const volScalarField& rhoField = forceSubM(0).rhoField();
    const volScalarField& nufField = forceSubM(0).nuField();

    scalar scaleDia3 = typeCG_[0]*typeCG_[0]*typeCG_[0];
    scalar cg = typeCG_[0];
    label partType = 1;

    //-------------Added for NN-------//
    if (nnDragCorrectorPtr_.valid()) {
        nnDragCorrectorPtr_().preLoop();
    }
    
    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            const label cellI = particleCloud_.cellIDs()[index][0];
            
            vector drag = vector::zero;
            vector dragExplicit = vector::zero;
            scalar dragCoefficient = 0;
            
            vector position = vector::zero;
            scalar voidfraction = 1.;
            vector Ufluid = vector::zero;

            if (cellI > -1) // particle Found
            {
                if(forceSubM(0).interpolation() )
                {
	            position     = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid       = UInterpolator_.interpolate(position,cellI);
                    if(voidfraction>1.00) voidfraction = 1.0;
                    if(voidfraction<0.10) voidfraction = 0.10;
                }
                else
                {
                    voidfraction = voidfraction_[cellI];
                    Ufluid       = U_[cellI];
                }

                if (multiTypes_)
                {
                    partType = particleCloud_.particleType(index);
                    if (partType > maxTypeCG_)
                    {
                        FatalError
                        << "Too few coarse-graining factors provided."
                        << abort(FatalError);
                    }
                    cg = typeCG_[partType - 1];
                    scaleDia3 = cg*cg*cg;
                }

                vector Us = particleCloud_.velocity(index);
                const scalar ds = 2*particleCloud_.radius(index);
                const scalar ds_scaled = ds/cg;
                const scalar rho = rhoField[cellI];
                const scalar nuf = nufField[cellI];

                const scalar Vs = ds*ds*ds*M_PI/6;
                const vector Ur = Ufluid-Us;
                const scalar magUr = mag(Ur);

                const scalar Rep = ds_scaled*voidfraction*magUr/(nuf+SMALL);
                const scalar CdMagUrLag = (24.0*nuf/(ds_scaled*voidfraction)) 
                                         *(1.0d+0.15*Foam::pow(Rep, 0.687));

                const scalar betaP 
                  = 0.75 * (
                              rho*voidfraction*CdMagUrLag
                              /
                              (ds_scaled*Foam::pow(voidfraction,2.65))
                           );
                
                scalar dragCoefficient = Vs*betaP;
                
                if (modelType_=="B")
                    dragCoefficient /= voidfraction;

                if (nnDragCorrectorPtr_.valid())
                {
                    const scalar Vcell = particleCloud_.mesh().V()[cellI];
                    const scalar dragCorrection =
                        nnDragCorrectorPtr_().calcDragCorrection(
                            cellI, voidfraction, position, Ur, Vcell, ds_scaled, nuf
                            );

                    dragCoefficient *= dragCorrection;
                }

                dragCoefficient *= scaleDia3*scaleDrag_;
                drag = dragCoefficient * Ur;
                dragExplicit = vector::zero;

                // explicitCorr
                for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                {
                    forceSubM(iFSub).explicitCorr(  drag, 
                                                    dragExplicit,
                                                    dragCoefficient,
                                                    Ufluid, U_[cellI], Us, UsField_[cellI],
                                                    forceSubM(iFSub).verbose()
                                                );
                }

                if(forceSubM(0).verbose() && index >=0 && index <2)
                {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "dprim = " << ds << endl;
		    Pout << "ds_scaled = " << ds_scaled << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "betaP = " << betaP << endl;
                    Pout << "drag = " << drag << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, drag);           //first entry must the be the force
                    vValues.setSize(vValues.size()+1, Ur);
                    sValues.setSize(sValues.size()+1, Rep); 
                    sValues.setSize(sValues.size()+1, betaP);
                    sValues.setSize(sValues.size()+1, voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
    }

    // write particle based data to global array
    forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);
    }// end loop particles
}
/*
double WenYuANNDrag::terminalVelocity(double voidfraction, double dp, double nuf, double rhof, double rhop, double g) const
{
    scalar u0 = dp*dp*fabs(rhof-rhop)*g/18./rhof/nuf;
    scalar Re = u0*dp/nuf;
    scalar res = 1.;
    scalar u = u0;
    scalar Fi(0);
    scalar CdSt(0);
    Info << "u0: " << u0 << endl;
    int i = 0;
    while ((res > 1.e-6) && (i<100))
    {
        Info << "Iteration " << i;
        u0 = u;
        Info << ", u0 = " << u0;
        CdSt = 24/Re;
        Info << ", CdSt = " << CdSt;
        Fi = F(voidfraction, Re);
        Info << ", F = ";
        u = sqrt(1.333333333*fabs(rhof-rhop)*g*dp
                  /(CdSt*voidfraction*Fi*rhof)
            );
        Info << ", u = " << u;
        Re = fabs(u)*dp/nuf*voidfraction;
        res = fabs((u-u0)/u);
        Info << "Res: " << res << endl;
        i++;
    }
    if (res > 1.e-6) {
        FatalError<< "Terminal velocity calculation diverged!" << endl;
    }

    return u;
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
