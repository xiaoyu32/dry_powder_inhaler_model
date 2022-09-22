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

#include "polarizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<PolarizationModel> PolarizationModel::New
(
    const dictionary& dict,
    cfdemCloud& sm
)
{
    Info<< "Selecting PolarizationModel "<<endl;
    //    << PolarizationModelType << endl;
    
    autoPtr<PolarizationModel> ptr;

    PolarizationModel* mm = new PolarizationModel( dict, sm );
    ptr.set( mm );
 
    Info<<"PolarizationModel selected."<<endl;

    return ptr;

    /*
    
    word PolarizationModelType
    (
        dict.lookup("PolarizationModel")
    );

    Info<< "Selecting PolarizationModel "
         << PolarizationModelType << endl;

    Info<<"Here we are!!!!!!!!!"<<endl;
    
    if( !dictionaryConstructorTablePtr_ )
    {
        Info<<"NULL dictionaryConstructorTablePtr! \n"<<endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(chargeDensityModelType);
    
    Info<<"Here we are!!!!!!!!!(2)"<<endl;
    
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "chargeDensityModelType::New(const dictionary&, const spray&) : "
            << endl
            << "    unknown chargeDensityModelType type "
            << chargeDensityModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid chargeDensityModelType types are :"
            << endl;
        Info<< dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }
    
    Info<<"ChargeDensityModel selected."<<endl;
    

    

    return autoPtr<ChargeDensityModel>(cstrIter()(dict,sm));
    */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
