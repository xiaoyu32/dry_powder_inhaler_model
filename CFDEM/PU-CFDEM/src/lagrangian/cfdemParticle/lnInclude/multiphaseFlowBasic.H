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


#ifndef multiphaseFlowBasic_H
#define multiphaseFlowBasic_H

#include "fvCFD.H"
#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiphaseFlowBasic Declaration
\*---------------------------------------------------------------------------*/

class multiphaseFlowBasic
{
public:
    bool    verbose_;

protected:

    // Protected data

public:

    // Constructor
    multiphaseFlowBasic();

    // Destructor
    ~multiphaseFlowBasic();
    
    //variables
    struct settlingParams
    {   
        double  u;              //the settling velocity of an isolated particle
        double  Re;             //the Re of an isolated particle
        double  Tref;           //u/g, the reference time scale
        double  Lref;           //u^2/g, the reference length scale
        double  Lchar;          //u^2/g*FrP^-0.5, a characteristic length scale
        double  Lchar2;          //u^2/g*FrP^-0.5, a characteristic length scale
        double  FrP;             //u^2/g/dp, the particle Froude number
    };
    
    settlingParams settling;
    
    // Member Functions
    double
    F_TangKuipers
              ( 
            double, 
            double
              );

    double
    F_WenYu
              ( 
            double, 
            double
              );
              
    double              
    F_KochHill
              ( 
            double, 
            double
              );

    void 
    uSettlingSphereSuspension( double dP, 
                               double rhoP, 
                               double g, 
                               double etaFluid, 
                               double rhoFluid, 
                               double phiP, 
                               int    dragLaw,
                               double &u,
                               double &Re
                               );
                               
    void 
    setupSettling(             double dP, 
                               double rhoP, 
                               double g, 
                               double etaFluid, 
                               double rhoFluid, 
                               int    dragLaw
                 );


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

