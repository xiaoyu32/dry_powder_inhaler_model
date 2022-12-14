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
    
Contributing Author
    Stefan Radl (Graz UT)
\*---------------------------------------------------------------------------*/


#include "multiphaseFlowBasicTK.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from scratch
multiphaseFlowBasicTK::multiphaseFlowBasicTK
(
)
:
    verbose_(false)
{
    settling.u        = 0.0;
    settling.Re       = 0.0;
    settling.Tref     = 0.0;
    settling.Lref     = 0.0;
    settling.Lchar    = 0.0;
    settling.Lchar2   = 0.0;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //
multiphaseFlowBasicTK::~multiphaseFlowBasicTK()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void 
multiphaseFlowBasicTK::setupSettling( double dP, 
                                    double rhoP, 
                                    double g, 
                                    double etaFluid, 
                                    double rhoFluid, 
                                    int    dragLaw
                                  )
{
    //Calculate the settling velocity of an isolated particle
    uSettlingSphereSuspension
                        ( 
                           dP, 
                           rhoP, 
                           g, 
                           etaFluid, 
                           rhoFluid, 
                           0, 
                           dragLaw,
                           settling.u,
                           settling.Re
                        );
    settling.Tref = settling.u / g;
    settling.Lref = settling.u*settling.u / g;
    settling.FrP  = settling.u*settling.u / g / dP;
    settling.Lchar = settling.Lref * pow(settling.FrP, -0.5);
    settling.Lchar2= settling.Lref * pow(settling.FrP, -0.666666666667);

    if(verbose_)
    {           
        Info << "\n *** Settling properties have been set up ***" << endl;
        Info << "  u:      " << settling.u << tab << "[m/s]" << endl;
        Info << "  Re:     " << settling.Re << endl;
        Info << "  Tref:   " << settling.Tref << tab << "[s]" << endl;
        Info << "  Lref:   " << settling.Lref << tab << "[m]" << endl;
        Info << "  FrP:    " << settling.FrP << endl;
        Info << "  Lchar:  " << settling.Lchar << tab << "[m]" <<  endl;
        Info << "  Lchar2: " << settling.Lchar2 << tab << "[m]" <<  endl << endl;
    }
    return;
}; //end function



void 
multiphaseFlowBasicTK::uSettlingSphereSuspension
                        ( 
                           double dP, 
                           double rhoP, 
                           double g, 
                           double etaFluid, 
                           double rhoFluid, 
                           double phiP, 
                           int    dragLaw,
                           double &u,
                           double &Re
                          )
{
//uSettlingSphere - Calculates the Settling Velocity of a sphere in a
//homogeneous suspension
//   uSettlingSphereSuspension( dP, rhoP, g, etaFluid, rhoFluid, phiP,
//   dragLaw )
// -------------------------------------
//   INPUT:
//   dP       particle diameter
//   rhoP     particle density
//   g        gravitational acceleration
//   etaFluid dyn. viscosity of fluid
//   rhoFluid density of the fluid
//   phiP     particle volume fraction
//   dragLaw  indicator for drag law
//               0Wen-Yu
//               1TangKuipers et al.
// -------------------------------------
//   OUTPUT:
//   u        settling velocity
//   Re       Reynolds number


scalar uInit  = dP * dP * (rhoP-rhoFluid) *g / (18 * etaFluid); //init with stokes settling velocity
scalar ReInit = fabs(uInit) * dP * rhoFluid / etaFluid;

scalar relError = 1.0;
scalar uOld(0.0);
uOld = uInit;
scalar cDStokes;
scalar F;

u     = uInit;
Re    = ReInit;

while(relError > 1.0e-6)
{
    
    // Calculate the drag coefficient
    cDStokes = 24/ Re; //Stokes drag coefficient
    if(dragLaw==0)      // Wen-Yu Drag law
        F = F_WenYu(    phiP, Re );
    else if(dragLaw==1)
        F = F_TangKuipers( phiP, Re );
    else if(dragLaw==2)
        F = F_KochHill( phiP, Re );
    else
        F = 1.;
    
    u    = sqrt( 4./3. * fabs(rhoP-rhoFluid) * g * dP 
            / ( cDStokes * (1.-phiP) * F * rhoFluid  ) 
             );
             
    Re  = (1.-phiP) * fabs(u) * dP * rhoFluid  / etaFluid;
    
    relError = fabs(uOld - u);
    uOld = u;

}

 return;

}; //end function


double
multiphaseFlowBasicTK::F_TangKuipers
          ( 
            double phiP, 
            double Re
          )
{
    double F;
    F    = 10. * phiP / (1.-phiP) / (1.-phiP) 
           + (1.-phiP)*(1.-phiP) * (1+1.5*sqrt(phiP)) 
	   + Re * ( 0.11 * phiP * ( 1. + phiP )
	           - 0.00456 / (1.-phiP) / (1.-phiP) / (1.-phiP) / (1.-phiP)
		   + ( 0.169 * (1.-phiP) + 0.0644 / (1.-phiP) / (1.-phiP) / (1.-phiP) / (1.-phiP) ) * 				     pow(Re,-0.343));

           /*+  0.413 * Re / 24. / sqr(1-phiP) 
            *( 1./(1.-phiP) + 3.0 * phiP * (1-phiP) + 8.4 * pow(Re,-0.343) )
            /( 1. + pow(10.,3.*phiP) * pow(Re,-(1.+4.*phiP)/2.) );*/

    return F;
}; //end function

double
multiphaseFlowBasicTK::F_WenYu
          ( 
            double phiP, 
            double Re
          )
{

    double F;

    F   = pow(1.-phiP,-3.65) 
        *( 1. + 0.15*pow(Re, 0.687) );


    return F;
}; //end function


double
multiphaseFlowBasicTK::F_KochHill
          ( 
            double phiP, 
            double Re
          )
{

    double F;
    double F0, F3;
    double voidfraction = 1. - phiP;
        
    //F0
    if(phiP < 0.4)
    {
              F0 = ( 1.+3.*sqrt((phiP)/2)+135./64.*phiP*log(phiP+1.e-16) 
                        +16.14*phiP 
                   ) / 
                   ( 1+0.681*phiP-8.48*phiP*phiP 
                              +8.16*phiP*phiP*phiP 
                   );
    }
    else
    {
              F0 = 10.*phiP/(voidfraction*voidfraction*voidfraction);
    }
        
    //F3
    F3 = 0.0673 + 0.212*phiP + 0.0232 / pow(voidfraction,5);
    
    F = voidfraction*(F0 + F3 * Re);

    return F;
}; //end function


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
 
