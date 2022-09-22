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

// Liquid volume in a bridge between particle i & js
double **liquidVol = new double*[maxNumberOfParticles];
for(int i = 0; i < maxNumberOfParticles; ++i) 
{
    liquidVol[i] = new double[3];
}

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

