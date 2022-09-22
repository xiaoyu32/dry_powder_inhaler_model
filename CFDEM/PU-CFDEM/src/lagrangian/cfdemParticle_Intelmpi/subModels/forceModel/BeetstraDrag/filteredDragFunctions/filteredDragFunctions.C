/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
                                Copyright 2013-     Graz University of Technology
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

// ***************************************************
// MODELS DESIGNED FOR THE BEETSTRA DRAG COEFFICIENT
// ***************************************************

inline double
BeetstraDrag::cCorrFunctionFilteredDrag(
                                          double k,     
                                          double aLimit,
                                          double aExponent,
                                          double phiP,
                                          double dParceldPrim
                                       ) const
{
    double phiPCrit = 0.55; //hard-coded particle critical volume fraction
                                      //for switching-off cCorrection
    double phiPStar = phiP / phiPCrit;
    if(phiPStar>1.0)  phiPStar = 1.0; 

    double aLimitStar = aLimit*phiPStar;
    if(aExponent>1.0) aLimitStar *=  pow(phiPStar,  aExponent-1);

    double result =    aLimitStar
                  +(1-aLimitStar)*exp(-k*(dParceldPrim-1));

    return  result;
}

inline double
BeetstraDrag::FfFilterFunc(
                            double LChar2,
                            double vCell,
                            double FroudeP
                          ) const
{
    //Function to calculate a dimensionless (inverse) filter size
    //that can be then used to correct the drag coefficient
    
    double result = 1.0;
    
    result = LChar2 / ( pow(vCell, 0.333333333333333) + SMALL);
    
    return result;
    
}

inline double
BeetstraDrag::fFuncFilteredDrag(
                   double    FfFilterPrime, 
                   double    localPhiP
                 ) const
{
    //Function to calculate a dimensionless (inverse) filter size
    //that can be then used to correct the drag coefficient
    
    //Constant within the function
    double   phiPCrit = 0.016f;
    double   critValueA = 90;
    double   breaks[6];
            breaks[0]= 0.00;
            breaks[1]= 0.100;
            breaks[2]= 0.180;
            breaks[3]= 0.250;
            breaks[4]= 0.400;
            breaks[5]= 0.550;
            
    double   coeff_x3[5];
    double   coeff_x2[5];
    double   coeff_x1[5];
    double   coeff_x0[5];
            coeff_x3[0] = -8049;
            coeff_x3[1] = 9431;
            coeff_x3[2] = -2523;
            coeff_x3[3] = 470;
            coeff_x3[4] = 0;
            
            coeff_x2[0] = 690.090;
            coeff_x2[1] = -1724.49;
            coeff_x2[2] = 538.880;
            coeff_x2[3] = 9.12;
            coeff_x2[4] = 0;
            
            coeff_x1[0] = 123.020;
            coeff_x1[1] = 19.5800;
            coeff_x1[2] = -75.2700;
            coeff_x1[3] = -36.9100;
            coeff_x1[4] = 0;
            
            coeff_x0[0] = 8.22;
            coeff_x0[1] = 19.3700;
            coeff_x0[2] = 14.7300;
            coeff_x0[3] = 11.2300;
            coeff_x0[4] = 7.49;


    //Calculation
    double  phiDiffSquared = 0.0f;
    double  phiDiffCubed    = 0.0f;
    double   a_Ff       = 1.0f;
    double   a_FfPrime  = 1.0f;
    int     rangeToUse = 4; //range to be used for large phiP values

    if(localPhiP > breaks[4])
    {
        rangeToUse = 4;
    }
    if(localPhiP < breaks[4])
    {
        rangeToUse = 3;
    }
    if(localPhiP < breaks[3])
    {
        rangeToUse = 2;
    }
    if(localPhiP < breaks[2])
    {
        rangeToUse = 1;
    }
    if(localPhiP < breaks[1])
    {
        rangeToUse = 0;
    }    
    if(localPhiP < phiPCrit)
    {
        rangeToUse = -1;
        a_Ff = critValueA;
    }
 
    if(rangeToUse>-1)
    {

         phiDiffSquared = (localPhiP-breaks[rangeToUse])
                                 *(localPhiP-breaks[rangeToUse]);
         phiDiffCubed   =  phiDiffSquared
                                 *(localPhiP-breaks[rangeToUse]);
        a_Ff = 
            phiDiffCubed                                  * coeff_x3[rangeToUse] 
          +phiDiffSquared                                * coeff_x2[rangeToUse] 
          +(localPhiP-breaks[rangeToUse])       * coeff_x1[rangeToUse] 
          +                                                      coeff_x0[rangeToUse];
    }
    
    a_FfPrime = 0.239102 * a_Ff; // 0.239 is the dp-to-LChar2 ratio of the system used for calibration
    
    double   result =  1.0 
                    / ( FfFilterPrime * a_FfPrime + sign(localPhiP - phiPCrit) + SMALL );
    
    return  result;
}           
                                       
inline double                                      
BeetstraDrag::hFuncFilteredDrag(
                   double   localPhiP
                 ) const
{
    //Constant within the function
    double maxPhiP   = 0.60;

    double   breaks[8];
            breaks[0]= 0.000;
            breaks[1]= 0.030;
            breaks[2]= 0.080;
            breaks[3]= 0.120;
            breaks[4]= 0.180;
            breaks[5]= 0.340;
            breaks[6]= 0.480;
            breaks[7]= 0.550;
            
    double   coeff_x3[8];
    double   coeff_x2[8];
    double   coeff_x1[8];
    double   coeff_x0[8];
            coeff_x3[0] = 0.0;  
            coeff_x3[1] = 253.63000;
            coeff_x3[2] = -789.56000;
            coeff_x3[3] = 310.8;
            coeff_x3[4] = 5.99;
            coeff_x3[5] = -9.12;
            coeff_x3[6] = -132.63000;
            coeff_x3[7] = 0;
            
            coeff_x2[0] =  0.0;
            coeff_x2[1] =  -4.41;
            coeff_x2[2] =  33.63;
            coeff_x2[3] = -61.11;
            coeff_x2[4] = -5.17;
            coeff_x2[5] = -2.29;
            coeff_x2[6] = -6.13;
            coeff_x2[7] = -225.2;
            
            coeff_x1[0] = 7.97;
            coeff_x1[1] = 4.64;
            coeff_x1[2] = 6.10;
            coeff_x1[3] = 5.01;
            coeff_x1[4] = 1.03;
            coeff_x1[5] = -0.170;
            coeff_x1[6] = -1.35;
            coeff_x1[7] = -2.3423;
            
            coeff_x0[0] =  0.0000;
            coeff_x0[1] =  0.239;
            coeff_x0[2] =  0.492;
            coeff_x0[3] =  0.739;
            coeff_x0[4] = 0.887;
            coeff_x0[5] = 0.943;
            coeff_x0[6] = 0.850;
            coeff_x0[7] = 0.6800;

    //Calculation
    double  phiDiffSquared = 0.0f;
    double  phiDiffCubed    = 0.0f;
    int     rangeToUse = 4; //range to be used for large phiP values
    if(localPhiP > maxPhiP )
    {
       rangeToUse = -1;
    }
    if(localPhiP > breaks[7])
    {
        rangeToUse = 7;
    }
    if(localPhiP < breaks[7])
    {
        rangeToUse = 6;
    }
    if(localPhiP < breaks[6])
    {
        rangeToUse = 5;
    }
    if(localPhiP < breaks[5])
    {
        rangeToUse = 4;
    }
    if(localPhiP < breaks[4])
    {
        rangeToUse = 3;
    }
    if(localPhiP < breaks[3])
    {
        rangeToUse = 2;
    }
    if(localPhiP < breaks[2])
    {
        rangeToUse = 1;
    }
    if(localPhiP < breaks[1])
    {
        rangeToUse = 0;
    }   

    double   result = 0.0;

    if(rangeToUse>-1)
    {
         phiDiffSquared = (localPhiP-breaks[rangeToUse])
                                 *(localPhiP-breaks[rangeToUse]);
         phiDiffCubed   =  phiDiffSquared
                                 *(localPhiP-breaks[rangeToUse]);

        result = 
            phiDiffCubed                                   * coeff_x3[rangeToUse] 
          +phiDiffSquared                                 * coeff_x2[rangeToUse] 
          +  (localPhiP-breaks[rangeToUse])      * coeff_x1[rangeToUse] 
          +                                                         coeff_x0[rangeToUse];
    }

    return result;
}           
