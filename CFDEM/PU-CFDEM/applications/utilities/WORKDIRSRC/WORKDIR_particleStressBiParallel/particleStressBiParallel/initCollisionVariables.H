// Collision parameters
int collisionModelI(0);

double* youngsModulus = new double[nParticleClass];
double* poissonsRatio = new double[nParticleClass];

double **coefficientRestitution = new double*[nParticleClass];
for(int i = 0; i < nParticleClass; ++i) 
{
    	coefficientRestitution[i] = new double[nParticleClass];
}	

double **coefficientFriction = new double*[nParticleClass];
for(int i = 0; i < nParticleClass; ++i) 
{
    	coefficientFriction[i] = new double[nParticleClass];
}

//double* coefficientRestitution = new double[nParticleClass*nParticleClass]; 
//double* coefficientFriction = new double[nParticleClass];

double k_n(0); // = new double[nParticleClass];
double k_t(0); // = new double[nParticleClass];
double gamma_n(0); // = new double[nParticleClass];
double gamma_t(0); // = new double[nParticleClass];
double mu_f(0); //= new double[nParticleClass];
double e_n(0); // = new double[nParticleClass];
double e_t(0); // = new double[nParticleClass];

double dt(0);

bool tangential_history=false;
if( calcCollision )
{
	if(ParticleStressPropsDict.found("tangential_history")) 	
	{
		tangential_history = true;
		Info << "Tangential history effect is active " << endl;
	}
}
bool liquid_transfer=false;
if(ParticleStressPropsDict.found("liquid_transfer")) liquid_transfer = true;	
double surf_tension(0);
double fluid_visc(0);

bool cohesion=false;
if(ParticleStressPropsDict.found("cohesion")) cohesion = true;	
double minimumDistanceVdW(0);
//double* cohEnergyDens = new double[nParticleClass];

double **cohEnergyDens = new double*[nParticleClass];
for(int i = 0; i < nParticleClass; ++i) 
{
    	cohEnergyDens[i] = new double[nParticleClass];
}

double *liq = new double[maxNumberOfParticles];

double **fcoll = new double*[maxNumberOfParticles];
for(int i = 0; i < maxNumberOfParticles; ++i) 
{
    fcoll[i] = new double[3];
}

double **ftan = new double*[maxNumberOfParticles];
for(int i = 0; i < maxNumberOfParticles; ++i) 
{
    ftan[i] = new double[3];
}

double **fcap = new double*[maxNumberOfParticles];
for(int i = 0; i < maxNumberOfParticles; ++i) 
{
    fcap[i] = new double[3];
}

double **fvisc = new double*[maxNumberOfParticles];
for(int i = 0; i < maxNumberOfParticles; ++i) 
{
    fvisc[i] = new double[3];
}

double **fcoh = new double*[maxNumberOfParticles];
for(int i = 0; i < maxNumberOfParticles; ++i) 
{
    fcoh[i] = new double[3];
}

// Liquid volume in a bridge between particle i & js
double **liquidVol = new double*[maxNumberOfParticles];
for(int i = 0; i < maxNumberOfParticles; ++i) 
{
    liquidVol[i] = new double[3];
}

/*
// Liquid bridge form
bool** first_touch = new bool*[maxNumberOfParticles]; 
for(int i = 0; i < maxNumberOfParticles; ++i) 
{
    first_touch[i] = new bool[maxNumberOfParticles];
}

for(int i = 0; i < maxNumberOfParticles; ++i) 
{
	for(int j = 0; j < maxNumberOfParticles; ++j) 
	{
		first_touch[i][j] = false;		
	}	
}

// Tangential displacement, particle-pair shear history
double*** delta_t = new double**[maxNumberOfParticles]; 
for(int i = 0; i < maxNumberOfParticles; ++i) 
{
    delta_t[i] = new double*[maxNumberOfParticles];
    for(int j = 0; j < maxNumberOfParticles; ++j) 
	{
        delta_t[i][j] = new double[3];
    }
}
*/

// Liquid bridge form
bool first_touch(false); 

// Tangential displacement, particle-pair shear history
double * delta_t = new double[3];
