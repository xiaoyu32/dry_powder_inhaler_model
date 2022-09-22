double **pos = new double*[MaxNumberOfParticles];
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    pos[i] = new double[3];
}	

double **vel = new double*[MaxNumberOfParticles];
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    vel[i] = new double[3];
}		

double **omega = new double*[MaxNumberOfParticles];
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    omega[i] = new double[3];
}	

double *liq = new double[MaxNumberOfParticles];
double *rad = new double[MaxNumberOfParticles];

double **fcoll = new double*[MaxNumberOfParticles];
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    fcoll[i] = new double[3];
}

double **ftan = new double*[MaxNumberOfParticles];
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    ftan[i] = new double[3];
}

double **fcap = new double*[MaxNumberOfParticles];
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    fcap[i] = new double[3];
}

double **fvisc = new double*[MaxNumberOfParticles];
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    fvisc[i] = new double[3];
}

double *density = new double[MaxNumberOfParticles];
double *mass = new double[MaxNumberOfParticles];

// Tangential displacement, particle-pair shear history
double*** delta_t = new double**[MaxNumberOfParticles]; 
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    delta_t[i] = new double*[MaxNumberOfParticles];
    for(int j = 0; j < MaxNumberOfParticles; ++j) 
	{
        delta_t[i][j] = new double[3];
    }
}

// Liquid volume in a bridge between particle i & js
double **liquidVol = new double*[MaxNumberOfParticles];
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    liquidVol[i] = new double[3];
}

// Liquid bridge form
bool** first_touch = new bool*[MaxNumberOfParticles]; 
for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
    first_touch[i] = new bool[MaxNumberOfParticles];
}

for(int i = 0; i < MaxNumberOfParticles; ++i) 
{
	for(int j = 0; j < MaxNumberOfParticles; ++j) 
	{
		first_touch[i][j] = false;		
	}	
}