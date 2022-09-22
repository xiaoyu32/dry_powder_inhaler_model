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

\*---------------------------------------------------------------------------*/

#include "CPCCellToCellStencil_mod.H"
#include "syncTools.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per point the neighbour data (= pointCells)
void Foam::CPCCellToCellStencil_mod::calcPointBoundaryData
(
    const boolList& isValidBFace,
    const labelList& boundaryPoints,
    Map<labelList>& neiGlobal
) const
{
    neiGlobal.resize(2*boundaryPoints.size());

    labelHashSet pointGlobals;

    forAll(boundaryPoints, i)
    {
        label pointI = boundaryPoints[i];

        neiGlobal.insert
        (
            pointI,
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointI],
                pointGlobals
            )
        );
    }

    syncTools::syncPointMap
    (
        mesh(),
        neiGlobal,
        unionEqOp(),
        Foam::dummyTransform()      // dummy transformation
    );
}


// Calculates per cell the neighbour data (= cell or boundary in global
// numbering). First element is always cell itself!
void Foam::CPCCellToCellStencil_mod::calcCellStencil
(
    labelListList& globalCellCells,
	int& filterwidth
) const
{
    // Calculate points on coupled patches
    labelList boundaryPoints(allCoupledFacesPatch()().meshPoints());


    // Mark boundary faces to be included in stencil (i.e. not coupled or empty)
    boolList isValidBFace;
    validBoundaryFaces(isValidBFace);


    // Swap pointCells for coupled points
    Map<labelList> neiGlobal;
    calcPointBoundaryData
    (
        isValidBFace,
        boundaryPoints,
        neiGlobal
    );

    globalCellCells.setSize(mesh().nCells());

    // Do coupled points first

    forAll(boundaryPoints, i)
    {
        label pointI = boundaryPoints[i];

        const labelList& pGlobals = neiGlobal[pointI];

        // Distribute to all pointCells
        const labelList& pCells = mesh().pointCells(pointI);

        forAll(pCells, j)
        {
            label cellI = pCells[j];

            // Insert pGlobals into globalCellCells
            merge
            (
                globalNumbering().toGlobal(cellI),
                pGlobals,
                globalCellCells[cellI]
            );
        }
    }
    neiGlobal.clear(); 

    // Do remaining points cells
    labelHashSet pointGlobals;
			
	int countIter = 0;
	
    for (label pointI = 0; pointI < mesh().nPoints(); pointI++)
    {
		labelList pGlobals
        (
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointI],
                pointGlobals
            )
        );
	

        const labelList& pCells = mesh().pointCells(pointI);

        forAll(pCells, j)
        {
            label cellI = pCells[j];
				
            merge
            (
                globalNumbering().toGlobal(cellI),
                pGlobals,
                globalCellCells[cellI]
            );
				//int iwidth = 1;
				//while ( iwidth < filterwidth/2 )
				for (int iwidth = 0; iwidth < 2; iwidth++)
				{
					const labelList& CelltoPoints = mesh().cellPoints()[cellI];
					
					forAll(CelltoPoints,pointII) 
					{
				            const labelList& pCells_sub = mesh().pointCells(CelltoPoints[pointII]);
				        
							forAll(pCells_sub, j)
					        {
					            label cellI_sub = pCells_sub[j];
			
					            merge
					            (
					                globalNumbering().toGlobal(cellI_sub),
					                pGlobals,
					                globalCellCells[cellI_sub]
					            );
		    
						    		// Info << " point I " << pointI << " " << " point II " << pointII << " " << mesh().nPoints() << endl;
								countIter++;				   						
		   						cellI = cellI_sub;		    		
							}
					 }
					 iwidth++;	 					
				}	
			
				    	
		}																										
 	
				   																	
    }	
	
	Info << " " << endl;
	Info << " Number of iteration " << countIter << endl;
    

    		
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::CPCCellToCellStencil_mod::CPCCellToCellStencil_mod
(	
	const polyMesh& mesh, 
	const int& FilterWidth
)
:
    cellToCellStencil(mesh)
{
    // Calculate per cell the (point) connected cells (in global numbering)
    labelListList globalCellCells;
    int filterwidth = FilterWidth;
	calcCellStencil(*this, filterwidth); 
}


// ************************************************************************* //



/*	
		//if ( pointList[pointI] < 0)
		//{
			
		labelList pGlobals
        (
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointI],
                pointGlobals
            )
        );
	
		//Info << " point " << pointI <<endl;

        const labelList& pCells = mesh().pointCells(pointI);

        forAll(pCells, j)
        {
            label cellI = pCells[j];
					
            merge
            (
                globalNumbering().toGlobal(cellI),
                pGlobals,
                globalCellCells[cellI]
            );
		
				for( label iwidth = 1; iwidth < filterwidth/2; iwidth++)
				{
				    for (label pointII = 0; pointII < mesh().nPoints(); pointII++)
				    {

				        const labelList& pCells = mesh().pointCells(pointII);

				        forAll(pCells, j)
				        {
				            label cellI = pCells[j];
					
				            merge
				            (
				                globalNumbering().toGlobal(cellI),
				                pGlobals,
				                globalCellCells[cellI]
				            );
				    
					    Info << " point I " << pointI << " " << " point II " << pointII << " " << mesh().nPoints() << endl;
								countIter++;
						}
				    }			
				}		
			
				    
	
		
		}		

*/


/*
		{			
			const labelList& CelltoPoints = mesh().cellPoints(cellI);
			
			forAll(CelltoPoints, CelltoPointsI)
			{
				const labelList& pCells2 = mesh().pointCells(CelltoPoints[CelltoPointsI]);
				
				forAll(pCells2, j)
			        {
			            label cellI2 = pCells2[j];
					
			            merge
			            (
			                globalNumbering().toGlobal(cellI),
			                pGlobals,
			                globalCellCells[cellI]
			            );
				   
				    merge
			            (
			                globalNumbering().toGlobal(cellI2),
			                pGlobals,
			                globalCellCells[cellI2]
			            );
				    
				    //Info << " point I " << pointI << " " << mesh().nPoints() << endl;							
				}				
			}		
		}
*/		
/*
	labelList pointList(mesh().nPoints());
    for (label pointI = 0; pointI < mesh().nPoints(); pointI++)
    {
		pointList[pointI] = -1;
	}		
	int countIter = 0;
	
    for (label pointI = 0; pointI < mesh().nPoints(); pointI++)
    {
        
		if ( pointList[pointI] < 0)
		{
			
		labelList pGlobals
        (
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointI],
                pointGlobals
            )
        );
	
		//Info << " point " << pointI <<endl;

        const labelList& pCells = mesh().pointCells(pointI);

        forAll(pCells, j)
        {
            label cellI = pCells[j];
					
            merge
            (
                globalNumbering().toGlobal(cellI),
                pGlobals,
                globalCellCells[cellI]
            );
		
				Info << " point " << pointI << " cellI " << cellI << endl;
				const labelList& CelltoPoints = mesh().cellPoints()[cellI];
						
				forAll(CelltoPoints,pointII) 
				{
					pointList[CelltoPoints[pointII]] = 0;
				}
		}				
					
				for( label iwidth = 1; iwidth < filterwidth/2; iwidth++)
				{
				    for (label pointII = 0; pointII < mesh().nPoints(); pointII++)
				    {

				        const labelList& pCells = mesh().pointCells(pointII);
					        forAll(pCells, j)
					        {
					            label cellI = pCells[j];
					
					            merge
					            (
					                globalNumbering().toGlobal(cellI),
					                pGlobals,
					                globalCellCells[cellI]
					            );
				    
						    		// Info << " point I " << pointI << " " << " point II " << pointII << " " << mesh().nPoints() << endl;
							
							}
					    }			
				}	
			
				    
	
		
		}														
 	
	countIter++;				   						
			    																	
    }	
	
	Info << " " << endl;
	Info << " Number of iteration " << countIter << endl;

*/