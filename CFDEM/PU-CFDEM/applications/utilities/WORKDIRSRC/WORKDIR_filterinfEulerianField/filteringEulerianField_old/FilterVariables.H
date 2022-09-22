// Minimum filter width 3X3X3
//const int minFilterWidth = 3;
// Maximum filter width 
//const int maxFilterWidth = 5;
// Filter width increment 
//const int FilterIncrement = 2;

//const int maxFilterWidth(readScalar(propsDict.lookup("maxFilterWidth")));	
//const int minFilterWidth(readScalar(propsDict.lookup("minFilterWidth")));	
//const int FilterIncrement(readScalar(propsDict.lookup("FilterIncrement")));

int maxFilterWidth(readScalar(propsDict.lookup("maxFilterWidth")));	
int minFilterWidth(readScalar(propsDict.lookup("minFilterWidth")));	
int FilterIncrement(readScalar(propsDict.lookup("FilterIncrement")));

Info << "" << endl;

if( minFilterWidth % 2 == 0 )
{
	FatalError<< "Even number for the smallest filter stencil width !!!" << abort(FatalError);	
}
Info << "The smallest filter width = " << minFilterWidth <<"X" << minFilterWidth << "X" << minFilterWidth << endl;

if( maxFilterWidth % 2 == 0 )
{
	FatalError<< "Even number for the largest filter stencil width !!!" << abort(FatalError);	
}
Info << "The largest  filter width = " << maxFilterWidth <<"X" << maxFilterWidth << "X" << maxFilterWidth << endl;

if( FilterIncrement % 2 != 0 )
{
	FatalError<< "Odd number for the filter width increment !!!" << abort(FatalError);	
}

if( minFilterWidth != maxFilterWidth )
{
	Info << "Filter width increment = " << FilterIncrement << endl;	
}

Info << "" << endl;

// To use less memory
maxFilterWidth   = maxFilterWidth/2-1;
minFilterWidth   = minFilterWidth/2-1;
FilterIncrement /= 2;
	
// Filter stencil list	    
labelListListList StencilListFilter(maxFilterWidth*mesh.cells().size());

// Global variables
volScalarField baralpf = voidfraction;
volVectorField tildeUf = U;
volScalarField barPg = p;

// Vector fields
vectorField DragForce(particleCloud.numberOfParticles());
vectorField ResolvedDragForce(particleCloud.numberOfParticles());

vectorField barBetai(particleCloud.numberOfParticles());
vectorField gii(particleCloud.numberOfParticles());

vectorField parcelbarBetai(particleCloud.numberOfParticles());
vectorField parcelgii(particleCloud.numberOfParticles());

// Stats  
scalarField                             filterScalar(maxFilterWidth+1,scalar(0.0));
Field <Field <scalar> >                 dimFilterScalar(3,filterScalar);

Field <Field <Field <scalar> > >        giiCondalppFilter(nBin,dimFilterScalar);
Field <Field <Field <scalar> > >        numbergiiCondalppFilter(nBin,dimFilterScalar);

int sizeNpartParticleCoarsening = npartParticleCoarsening.size();

Field <Field <Field <Field <scalar> > > >   NparcelgiiCondalppFilter(sizeNpartParticleCoarsening+1,giiCondalppFilter);
Field <Field <Field <Field <scalar> > > >   NparcelnumbergiiCondalppFilter(sizeNpartParticleCoarsening+1,numbergiiCondalppFilter);

//Field <Field <Field <scalar> > >    NparcelgiiCondalppFilter(nBin,dimFilterScalar);
//Field <Field <Field <scalar> > >    NparcelnumbergiiCondalppFilter(nBin,dimFilterScalar);

// Initialisation				
for(int FilterWidth = minFilterWidth; FilterWidth <= maxFilterWidth; FilterWidth += FilterIncrement )
{
	for(int i = 0; i < nBin; i++ )
	{
		for(int j=0;j<3;j++) 
		{				  
 				      giiCondalppFilter[i][j][FilterWidth] = 0;
				numbergiiCondalppFilter[i][j][FilterWidth] = 0; 
				
			//NparcelgiiCondalppFilter[i][j][FilterWidth] = 0;
			//NparcelnumbergiiCondalppFilter[i][j][FilterWidth] = 0;				
				
				for(int k=0;k<=sizeNpartParticleCoarsening;k++)
				{
					NparcelgiiCondalppFilter[k][i][j][FilterWidth] = 0;
					NparcelnumbergiiCondalppFilter[k][i][j][FilterWidth] = 0;
				}
		}
											
	} 
}



