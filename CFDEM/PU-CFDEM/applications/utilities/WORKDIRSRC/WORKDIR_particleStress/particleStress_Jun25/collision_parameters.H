// Collision parameters 
double youngsModulus(0);
double poissonsRatio(0);
double coefficientRestitution(0); 
double coefficientFriction(0);
double k_n(0);
double k_t(0);
double gamma_n(0);
double gamma_t(0);
double mu_f(0);
double e_n(0);
double e_t(0);

double dt(0);

bool tangential_history=false;
if(ParticleStressPropsDict.found("tangential_history")) 	
{
	tangential_history = true;
	Info << "Tangential history effect is active " << endl;
}
bool liquid_transfer=false;
if(ParticleStressPropsDict.found("liquid_transfer")) liquid_transfer = true;	
double surf_tension(0);
double fluid_visc(0);

bool cohesion=false;
if(ParticleStressPropsDict.found("cohesion")) cohesion = true;	
double minimumDistanceVdW(0);
double cohEnergyDens(0);

fileName DEMinputFile(mesh.time().path()/DEMinputFileRelativePath/DEMinputFilename);
// Set file pointer
string DEMinputPointer=string(DEMinputFile);
const char * DEMFilePath=DEMinputPointer.c_str();
ifstream inputPtr(DEMFilePath);

// Collision model
int collisionModelI(0); // Hertz (default)
//word collisionModel("hooke");
const char * search="hooke";
string line;
while(getline(inputPtr, line))
{
	if (line.find(search, 0) != string::npos)
	{
		collisionModelI = 1;	
	}
}
inputPtr.close();

// Re-open the file
inputPtr.open(DEMFilePath);

if(collisionModelI==0)
{
	word name("youngsModulus");
	string just_read = " ";
	while(just_read.compare(name) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;			// Skip the phrase 
	inputPtr >> youngsModulus;
	word name2("poissonsRatio");
	while(just_read.compare(name2) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;                 // Skip the phrase 
	inputPtr >> poissonsRatio;
	word name3("coefficientRestitution");
	while(just_read.compare(name3) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;                 // Skip the phrase 
	inputPtr >> just_read;                 // Skip the phrase
	inputPtr >> coefficientRestitution;
	word name4("coefficientFriction");
	while(just_read.compare(name4) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;                 // Skip the phrase
	inputPtr >> just_read;                 // Skip the phrase 
	inputPtr >> coefficientFriction;

	Info << " " << endl;
	Info << "Collision model: Hertz" << endl;
	Info << "Collision parameters" << endl;
	Info << "Young's modulus    = " << youngsModulus << endl;        
	Info << "Poisson's ratio    = " << poissonsRatio << endl;                                   
	Info << "Restitution coeff. = " << coefficientRestitution << endl;
	Info << "Friction coeff.    = " << coefficientFriction << endl;
	
	if(cohesion)
	{
		word name5("minimumDistanceVdW");
		while(just_read.compare(name5) != 0)  inputPtr >> just_read;
		inputPtr >> just_read;                 // Skip the phrase
		inputPtr >> minimumDistanceVdW;		
		word name6("cohesionEnergyDensity");
		while(just_read.compare(name6) != 0)  inputPtr >> just_read;
		inputPtr >> just_read;                 // Skip the phrase
		inputPtr >> just_read;                 // Skip the phrase
		inputPtr >> cohEnergyDens;		
	
		Info << " " << endl;
		Info << "vdW Force is active   " << endl;                                   
		Info << "Min. distance       = " << minimumDistanceVdW << endl;
		Info << "Coh. energy density = " << cohEnergyDens << endl;

	}	

}

else if (collisionModelI==1)
{
	word name("kn");
	string just_read = " ";
	while(just_read.compare(name) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;			// Skip the phrase 
	inputPtr >> just_read;                 // Skip the phrase
	inputPtr >> k_n;
	word name2("kt");
	while(just_read.compare(name2) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;                 // Skip the phrase 
	inputPtr >> just_read;                 // Skip the phrase
	inputPtr >> k_t;
	word name3("gamman");
	while(just_read.compare(name3) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;                 // Skip the phrase 
	inputPtr >> just_read;                 // Skip the phrase
	inputPtr >> gamma_n;
	word name4("gammat");
	while(just_read.compare(name4) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;                 // Skip the phrase
	inputPtr >> just_read;                 // Skip the phrase 
	inputPtr >> gamma_t;

	Info << " " << endl;
	Info << "Collision model: Hooke" << endl;
	Info << "Collision parameters" << endl;
	Info << "Normal spring constant      = " << k_n << endl;        
	Info << "Tangential spring constant  = " << k_t << endl;                                   
	Info << "Normal damping constant     = " << gamma_n << endl;
	Info << "Tangential damping constant = " << gamma_t << endl;

	if(cohesion)
	{
		word name5("minimumDistanceVdW");
		while(just_read.compare(name5) != 0)  inputPtr >> just_read;
		inputPtr >> just_read;                 // Skip the phrase
		inputPtr >> minimumDistanceVdW;		
		word name6("cohesionEnergyDensity");
		while(just_read.compare(name6) != 0)  inputPtr >> just_read;
		inputPtr >> just_read;                 // Skip the phrase
		inputPtr >> just_read;                 // Skip the phrase
		inputPtr >> minimumDistanceVdW;	
		
		Info << " " << endl;
		Info << "vdW Force is active   " << endl;                                   
		Info << "Min. distance       = " << minimumDistanceVdW << endl;
		Info << "Coh. energy density = " << cohEnergyDens << endl;			
	}

}
else
{
	FatalError << " Define collision model in ParticleStressProbs " << abort(FatalError);
}	

// Define & Init collision variables
#include "init_collision_variables.H" 


