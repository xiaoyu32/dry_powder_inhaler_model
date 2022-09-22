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

// Eulerian Particle Velocity (Initialize with fluid velocity)
//volVectorField Up = U;
volVectorField Up
(
	IOobject
	(
		"Up",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionSet(0,1,-1,0,0,0,0)
);	
	
volVectorField tildeUp = U;

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
scalarField                             	filterScalar(maxFilterWidth+1,scalar(0.0));

Field <Field <scalar> >                 	dimFilterScalar(3,filterScalar);

Field <Field <scalar> >                 	taupVelSlipCondalppFilter(nBin,filterScalar);
Field <Field <scalar> >                 	numberCondalppFilter(nBin,filterScalar);

Field <Field <Field <scalar> > > 		NparcelTaupVelSlipCondalppFilter(npartParticleCoarsening.size(),taupVelSlipCondalppFilter);
Field <Field <Field <scalar> > >    		NparcelnumberCondalppFilter(npartParticleCoarsening.size(),numberCondalppFilter);

Field <Field <Field <scalar> > >        	giiCondalppFilter(nBin,dimFilterScalar);
Field <Field <Field <scalar> > >        	velSlipCondalppFilter(nBin,dimFilterScalar);

Field <Field <Field <Field <scalar> > > >   	NparcelgiiCondalppFilter(npartParticleCoarsening.size(),giiCondalppFilter);
Field <Field <Field <Field <scalar> > > >   	NparcelVelSlipCondalppFilter(npartParticleCoarsening.size(),velSlipCondalppFilter);

// Fluid coarsening joint gii[z]
Field <Field <Field <scalar> > >        	VelSlipJointgiiCondalppFilter(nBin,numberCondalppFilter);
Field <Field <Field <scalar> > > 		TaupVelSlipJointgiiCondalppFilter(nBin,numberCondalppFilter);
Field <Field <Field <scalar> > > 		numberVelSlipJointgiiCondalppFilter(nBin,numberCondalppFilter);
Field <Field <Field <scalar> > > 		numberTaupVelSlipJointgiiCondalppFilter(nBin,numberCondalppFilter);

// Parcel coarsening joint gii[z]
Field <Field <Field <Field <scalar> > > >      	NparcelVelSlipJointgiiCondalppFilter(npartParticleCoarsening.size(),VelSlipJointgiiCondalppFilter);
Field <Field <Field <Field <scalar> > > >  	NparcelTaupVelSlipJointgiiCondalppFilter(npartParticleCoarsening.size(),TaupVelSlipJointgiiCondalppFilter);
Field <Field <Field <Field <scalar> > > >   	NparcelnumberVelSlipJointgiiCondalppFilter(npartParticleCoarsening.size(),numberVelSlipJointgiiCondalppFilter);
Field <Field <Field <Field <scalar> > > >   	NparcelnumberTaupVelSlipJointgiiCondalppFilter(npartParticleCoarsening.size(),numberTaupVelSlipJointgiiCondalppFilter);

