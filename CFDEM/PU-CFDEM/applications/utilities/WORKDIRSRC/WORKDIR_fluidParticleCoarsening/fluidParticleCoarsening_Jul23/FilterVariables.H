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
	dimensionedVector( "zero", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0) )
);	
	
volVectorField tildeUp 
(
	IOobject
	(
		"tildeUp",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	Up
);

volVectorField tildeUf 
(
	IOobject
	(
		"tildeUf",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	U
);

// Global variables
volScalarField baralpf
(
	IOobject
	(
		"baralpf",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	voidfraction
);


volScalarField barPg
(
	IOobject
	(
		"barPg",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	p
);

// Vector fields
vectorField DragForce(particleCloud.numberOfParticles());
vectorField ResolvedDragForce(particleCloud.numberOfParticles());

vectorField barBetai(particleCloud.numberOfParticles());
vectorField gii(particleCloud.numberOfParticles());

vectorField parcelbarBetai(particleCloud.numberOfParticles());
vectorField parcelgii(particleCloud.numberOfParticles());

// Mapped Drag force
volVectorField MappedDragForce
(
    IOobject
    (
        "MappedDragForce",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
); 

// Eulerian drag force 
volVectorField EulerianDragForce
(
    IOobject
    (
        "EulerianDragForce",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);
// Eulerian resolved drag force
volVectorField ResolvedEulerianDragForce
(
    IOobject
    (
        "ResolvedEulerianDragForce",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

// Eulerian correction factor in vector array
vectorField Eulerian_gii(mesh.nCells());

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

// Eulerian filtering
Field <Field <Field <scalar> > >        	Eulerian_giiCondalppFilter(nBin,dimFilterScalar);
Field <Field <Field <scalar> > >        	Eulerian_velSlipCondalppFilter(nBin,dimFilterScalar);

Field <Field <scalar> >                 	Eulerian_taupVelSlipCondalppFilter(nBin,filterScalar);
Field <Field <scalar> >                 	Eulerian_numberCondalppFilter(nBin,filterScalar);

Field <Field <Field <scalar> > >        	Eulerian_VelSlipJointgiiCondalppFilter(nBin,numberCondalppFilter);
Field <Field <Field <scalar> > > 		Eulerian_TaupVelSlipJointgiiCondalppFilter(nBin,numberCondalppFilter);
Field <Field <Field <scalar> > > 		Eulerian_numberVelSlipJointgiiCondalppFilter(nBin,numberCondalppFilter);
Field <Field <Field <scalar> > > 		Eulerian_numberTaupVelSlipJointgiiCondalppFilter(nBin,numberCondalppFilter);
