
// Collision model
collisionModelI = 0; // Hertz (default)	

//if(Pstream::master()) //Write only if master
//{
//
const dictionary twoWayMPIPropsDict(dict.subDict("twoWayMPIProps"));	
// DEM input file parameters
fileName DEMinputFilename(twoWayMPIPropsDict.lookup("liggghtsPath"));

// Check DEM input file exist or not
//if( !isFile(mesh.time().path()/DEMinputFilename) )
if( !isFile(mesh.time()/DEMinputFilename) )
{
	FatalError << " DEM input file does not exist " << abort(FatalError);
}
else
{
	Info << " " << endl;
	Info << "DEM input file is " << mesh.time().path()/DEMinputFilename << endl;
}

fileName DEMinputFile(mesh.time().path()/DEMinputFilename);
// Set file pointer
string DEMinputPointer=string(DEMinputFile);
const char * DEMFilePath=DEMinputPointer.c_str();
ifstream inputPtr(DEMFilePath);

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
	Info << " " << endl;
	Info << "Collision model: Hertz" << endl;

	int iPartClass = 0;	// Mono-disperse 

	word name("youngsModulus");
	string just_read = " ";
	while(just_read.compare(name) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;			// Skip the phrase 
	//inputPtr >> youngsModulus;
	inputPtr >> youngsModulus[iPartClass];
	if( bidisperse ) inputPtr >> youngsModulus[iPartClass+1];

	word name2("poissonsRatio");
	while(just_read.compare(name2) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;                 // Skip the phrase 
	//inputPtr >> poissonsRatio;
	inputPtr >> poissonsRatio[iPartClass];
	if( bidisperse ) inputPtr >> poissonsRatio[iPartClass+1];	

	word name3("coefficientRestitution");
	while(just_read.compare(name3) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;                 // Skip the phrase 
	inputPtr >> just_read;                 // Skip the phrase
	//inputPtr >> coefficientRestitution;
	inputPtr >> coefficientRestitution[iPartClass][iPartClass];
	if( bidisperse ) 
	{
		inputPtr >> coefficientRestitution[iPartClass][iPartClass+1];
		inputPtr >> coefficientRestitution[iPartClass+1][iPartClass];
		inputPtr >> coefficientRestitution[iPartClass+1][iPartClass+1];		
	}

	word name4("coefficientFriction");
	while(just_read.compare(name4) != 0)  inputPtr >> just_read;
	inputPtr >> just_read;                 // Skip the phrase
	inputPtr >> just_read;                 // Skip the phrase 
	//inputPtr >> coefficientFriction;
	inputPtr >> coefficientFriction[iPartClass][iPartClass];
	if( bidisperse ) 
	{
		inputPtr >> coefficientFriction[iPartClass][iPartClass+1];
		inputPtr >> coefficientFriction[iPartClass+1][iPartClass];
		inputPtr >> coefficientFriction[iPartClass+1][iPartClass+1];	
	}

	if( !bidisperse ) 
	{
		Info << " " << endl;
		Info << "Particle class = " << iPartClass << " collision parameters" << endl;
		Info << "Young's modulus    = " << youngsModulus[iPartClass] << endl;        
		Info << "Poisson's ratio    = " << poissonsRatio[iPartClass] << endl;                                   
		Info << "Restitution coeff. = " << coefficientRestitution[iPartClass][iPartClass] << endl;
		Info << "Friction coeff.    = " << coefficientFriction[iPartClass][iPartClass] << endl;
	}	
	else
	{
		iPartClass = 0;
		Info << " " << endl;
		Info << "Particle class = " << iPartClass << " collision parameters" << endl;
		Info << "Young's modulus    = " << youngsModulus[iPartClass] << endl;        
		Info << "Poisson's ratio    = " << poissonsRatio[iPartClass] << endl;                                   
		Info << "Restitution coeff. = " << coefficientRestitution[iPartClass][iPartClass] << " " << coefficientRestitution[iPartClass][iPartClass+1] << endl;
		Info << "Friction coeff.    = " << coefficientFriction[iPartClass][iPartClass]    << " " << coefficientFriction[iPartClass][iPartClass+1] << endl;		

		iPartClass = 1;
		Info << " " << endl;
		Info << "Particle class = " << iPartClass << " collision parameters" << endl;
		Info << "Young's modulus    = " << youngsModulus[iPartClass] << endl;        
		Info << "Poisson's ratio    = " << poissonsRatio[iPartClass] << endl;                                   
		Info << "Restitution coeff. = " << coefficientRestitution[iPartClass][iPartClass] << " " << coefficientRestitution[iPartClass][iPartClass-1] << endl;
		Info << "Friction coeff.    = " << coefficientFriction[iPartClass][iPartClass]    << " " << coefficientFriction[iPartClass][iPartClass-1] << endl;					

	}	

	if(cohesion)
	{
		// Modify later minDist & smin for updated version of LIGGHTS

		word name5("minimumDistanceVdW");
		while(just_read.compare(name5) != 0)  inputPtr >> just_read;
		inputPtr >> just_read;                 // Skip the phrase
		inputPtr >> minimumDistanceVdW;	

		word name6("cohesionEnergyDensity");
		while(just_read.compare(name6) != 0)  inputPtr >> just_read;
		inputPtr >> just_read;                 // Skip the phrase
		inputPtr >> just_read;                 // Skip the phrase
		//inputPtr >> cohEnergyDens;		
		inputPtr >> cohEnergyDens[iPartClass][iPartClass];
		if( bidisperse ) 
		{
			inputPtr >> cohEnergyDens[iPartClass][iPartClass+1];
			inputPtr >> cohEnergyDens[iPartClass+1][iPartClass];
			inputPtr >> cohEnergyDens[iPartClass+1][iPartClass+1];				
		}

		Info << " " << endl;
		Info << "vdW Force is active   " << endl;                                   
		Info << "Min. distance       = " << minimumDistanceVdW << endl;

		if( !bidisperse)
		{
			Info << "Coh. energy density = " << cohEnergyDens[iPartClass][iPartClass+1] << endl;
		}
		else
		{
			iPartClass = 0;
			Info << "Particle class = " << iPartClass << " cohesion energy " << cohEnergyDens[iPartClass][iPartClass] << " " << cohEnergyDens[iPartClass][iPartClass+1] << endl;
			iPartClass = 1;
			Info << "Particle class = " << iPartClass << " cohesion energy " << cohEnergyDens[iPartClass][iPartClass] << " " << cohEnergyDens[iPartClass][iPartClass-1] << endl;
		}
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

//}
// Define & Init collision variables
//#include "init_collision_variables.H" 


