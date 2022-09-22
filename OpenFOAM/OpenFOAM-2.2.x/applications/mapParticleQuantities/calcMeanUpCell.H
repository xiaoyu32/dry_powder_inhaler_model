// Number of particles in each cell
volScalarField nPCell
(
     IOobject
     (
        "nP",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,0,0,0), scalar(0) )
);
   
for(int ii = 0; ii < nP[me]; ii++)
{							
        label cellI = cellIDs[ii];
        if( cellI>-1 ) 
	{
		nPCell[cellI] += 1.0;
		meanVelCell[cellI] += velocities[me][ii];	
	}
}
// Normalize velocity at each cell
forAll(mesh.cells(),cellI)
{
    if(nPCell[cellI]!=0)
    {
        for(int idir=0; idir<3; idir++) meanVelCell[cellI][idir]/= nPCell[cellI];
    }else
    {
        for(int idir=0; idir<3; idir++) meanVelCell[cellI][idir] = 0.;
    }
}

// Write into files
nPCell.write();
meanVelCell.write();
