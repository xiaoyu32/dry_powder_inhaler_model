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

#include "oneWayDump.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneWayDump, 0);

addToRunTimeSelectionTable
(
    dataExchangeModel,
    oneWayDump,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
oneWayDump::oneWayDump
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dataExchangeModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    filename_(propsDict_.lookup("couplingFilename")),
    relativePath_(propsDict_.lookup("relativePath"))
{
    readDEMtsfromDict(propsDict_);

    // set max nr of particles from dict
    maxNumberOfParticles_ = readScalar(propsDict_.lookup("maxNumberOfParticles"));
    setNumberOfParticles(maxNumberOfParticles_);

    // make a const char* from word
    string HH=string(filename_);
    charFilename_=HH.c_str();

    Info << "relativePath_" << relativePath_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oneWayDump::~oneWayDump()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void oneWayDump::getData
(
    word name,
    word type,
    double ** const& field,
    label step
) const
{

    // get path to particle VTK files
    char index[100];
    sprintf(index, charFilename_, step);
    //fileName H(particleCloud_.mesh().time().path()/".."/"DEM"/"post"/index);
    //fileName H(particleCloud_.mesh().time().path()/relativePath_/index);
    //_AO_09/19/2014 - Able to run parallel post-processing tools
        fileName H
        (
           particleCloud_.mesh().time().processorCase()
         ? particleCloud_.mesh().time().path()/".."/relativePath_/index
         : particleCloud_.mesh().time().path()/relativePath_/index
        );
    //_AO_09/19/2014 - Able to run parallel post-processing tools
    Info << "opening file: " <<H << endl;

    // set file pointer
    string HH=string(H);
    const char * paricleFilePath=HH.c_str();
    ifstream* inputPtr;
    inputPtr = new ifstream(paricleFilePath);
    
    // Read data
    string just_read = " ";
    int just_read_int;
    scalar just_read_scalar;
    label nP;
    
    double ** fieldDummy;
    particleCloud_.dataExchangeM().allocateArray(fieldDummy,0.,16);
    
    for(label ii=0; ii<2; ii++)
    {
    	*inputPtr >> just_read;   // ITEM: TIMESTEP
    }    
    *inputPtr >> just_read_int;   // number    
    for(label ii=0; ii<4; ii++)
    {
    	*inputPtr >> just_read;   // ITEM: NUMBER OF ATOMS
    }    
    *inputPtr >> nP;   // number
    
     // give nr of particles to cloud
     setNumberOfParticles(nP);

     // re-allocate arrays of cloud
     particleCloud_.reAllocArrays();
	    
    for(label ii=0; ii<6; ii++)
    {
    	*inputPtr >> just_read;   // ITEM: BOX BOUNDS pp pp pp
    }
    for(label ii=0; ii<2; ii++)
    {
    	*inputPtr >> just_read_scalar;   // number
    }
    for(label ii=0; ii<2; ii++)
    {
    	*inputPtr >> just_read_scalar;   // number
    } 
    for(label ii=0; ii<2; ii++)
    {
    	*inputPtr >> just_read_scalar;   // number
    }           
    for(label ii=0; ii<17; ii++)
    {
    	*inputPtr >> just_read;   // ITEM: ATOMS id type x y z vx vy vz fx fy fz omegax omegay omegaz radius 
    }
    
    for(int index = 0; index < nP; ++index)
    {
        *inputPtr >> fieldDummy[index][0]
	          >> fieldDummy[index][1]
		  >> fieldDummy[index][2]
	          >> fieldDummy[index][3]
		  >> fieldDummy[index][4]
	          >> fieldDummy[index][5]
		  >> fieldDummy[index][6]		  
	          >> fieldDummy[index][7]
		  >> fieldDummy[index][8]
	          >> fieldDummy[index][9]
		  >> fieldDummy[index][10]
	          >> fieldDummy[index][11]
		  >> fieldDummy[index][12]
		  >> fieldDummy[index][13]
	          >> fieldDummy[index][14];
		  
/*
        Info << fieldDummy[index][0] << " "
	          << fieldDummy[index][1] << " "
		  << fieldDummy[index][2] << " "
	          << fieldDummy[index][3] << " "
		  << fieldDummy[index][4] << " "
	          << fieldDummy[index][5] << " "
		  << fieldDummy[index][6] << " "		  
	          << fieldDummy[index][7] << " "
		  << fieldDummy[index][8] << " "
	          << fieldDummy[index][9] << " "
		  << fieldDummy[index][10] << " "
	          << fieldDummy[index][11] << " "
		  << fieldDummy[index][12] << " "
		  << fieldDummy[index][13] << " "
	          << fieldDummy[index][14] <<nl;	  
*/
    }    


    if (type == "scalar-atom")
    {

	if (name == "type")
	{
            for(int index = 0; index < nP; ++index)
    	    {
	       field[index][0] = fieldDummy[index][1];
            }
	}
	else if (name == "radius" )
	{
            for(int index = 0; index < nP; ++index)
    	    {
	       field[index][0] = fieldDummy[index][14];
            }
	}

	// clean up inputStream
	delete inputPtr;
    } else if (type == "vector-atom")
    {

	if (name == "x")
	{
            for(int index = 0; index < nP; ++index)
    	    {
	       field[index][0] = fieldDummy[index][2];
	       field[index][1] = fieldDummy[index][3];
	       field[index][2] = fieldDummy[index][4];
            }        
	}
	else if(name == "v")
	{
            for(int index = 0; index < nP; ++index)
    	    {
	       field[index][0] = fieldDummy[index][5];
	       field[index][1] = fieldDummy[index][6];
	       field[index][2] = fieldDummy[index][7];
	    }
	}
	else if(name == "force")
	{
            for(int index = 0; index < nP; ++index)
    	    {
	       field[index][0] = fieldDummy[index][8];
	       field[index][1] = fieldDummy[index][9];
	       field[index][2] = fieldDummy[index][10];
	    }
	}    
	else if(name == "omega")
	{
            for(int index = 0; index < nP; ++index)
    	    {
	       field[index][0] = fieldDummy[index][11];
	       field[index][1] = fieldDummy[index][12];
	       field[index][2] = fieldDummy[index][13];
	    }
	}
	// clean up inputStream
	delete inputPtr;
    }
    else
    {
	Info << "unknown type in getData!!!" << endl;
    }
    

}

void oneWayDump::giveData
(
    word name,
    word type,
    double ** const& field,
    const char* datatype
) const
{
    // do nothing
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
