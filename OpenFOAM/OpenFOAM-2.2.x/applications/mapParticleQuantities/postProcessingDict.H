// Post-processing dictionary 
IOdictionary couplingProperties
(
    IOobject
    (
        "couplingPropertiesMapParticleQuantities",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Define couplingProperties dictionary
const dictionary dict(couplingProperties);

// Define postprocessing dictionary
const dictionary postProcessingDict(dict.subDict("mapParticleQuantitiesProps"));

// Filename for dump files
word filename(postProcessingDict.lookup("couplingFilename"));
fileName relativePath(postProcessingDict.lookup("relativePath"));


word outputFilename(postProcessingDict.lookup("outputFilename"));

// DEM timestep
const scalar DEMts(readScalar(postProcessingDict.lookup("DEMts")));

// Debugging example particles
bool verbose(false);
labelList exList(1,1);
if(postProcessingDict.found("verbose"))
{
        verbose = true;
        exList = labelList(postProcessingDict.lookup("exList"));
}

// Parallel computation
int nProcs(1);
if(Pstream::parRun()) nProcs = Pstream::nProcs();

// OpenFOAM searching algorithm
Switch treeSearch(postProcessingDict.lookup("treeSearch"));
//meshSearch searchEngine(U.mesh(),polyMesh::FACEPLANES);
meshSearch searchEngine(mesh,polyMesh::FACEPLANES);

// Interpolation of Us OR Us vs Us*
bool interpolation(false);
if(postProcessingDict.found("interpolation")) 
{
	interpolation = true;
	Info << "Us will be interpolated" << endl;
}
bool useUsstar(false);
if(postProcessingDict.found("useUsstar")) 
{
	useUsstar = true;
	Info << "meanUp will be used instead of Us" << endl;	
}	

// Max. number of cells per particle
int maxCellPerPart(30);
if (postProcessingDict.found("maxCellPerPart"))
{
   maxCellPerPart = readScalar(postProcessingDict.lookup("maxCellPerPart"));
   Info << tab << "Max. number of cells per particle " << maxCellPerPart
             << " for mapping & interpolation" << endl;
}

// Min. alpha for mapping
scalar alphaMin(0.3); 

// Not weighting
bool notWeighthing(false);
if(postProcessingDict.found("notWeighting")) 
{
	notWeighthing = true;
	Info << "\nWeigthing is not used" << endl;
	maxCellPerPart = 1;
}


