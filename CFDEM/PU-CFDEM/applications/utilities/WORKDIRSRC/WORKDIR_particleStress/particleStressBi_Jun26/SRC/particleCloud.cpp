#include "particleCloud.h"

void particleCloud::set_values (int x, int y) 
{
  width = x;
  height = y;
}

int particleCloud::setNumberOfParticles (int ii) 
{
  nP = ii;
  return nP;
}

double* particleCloud::position(int index)
{
    double *pos = new double[3];
    for(int i=0;i<3;i++) pos[i] = positions[index][i];
    return pos;
}

double* particleCloud::velocity(int index)
{
    double *vel = new double[3];
    for(int i=0;i<3;i++) vel[i] = velocities[index][i];
    return vel;
}

double* particleCloud::omega(int index)
{
    double *ome = new double[3];
    for(int i=0;i<3;i++) ome[i] = omegas[index][i];
    return ome;
}

double particleCloud::liquidOn(int index) 
{
  	double liq;
  	liq = liquidOns[index];
  	return liq; 
}

double particleCloud::radius(int index) 
{
  	double rad;
  	rad = rads[index];
  	return rad; 
}

// Collisional force
double* particleCloud::force(int index) 
{
    double *force = new double[3];
    for(int i=0;i<3;i++) force[i] = forces[index][i];
    return force;
}

// Tangential collisional force
double* particleCloud::forceTan(int index) 
{
    double *forceT = new double[3];
    for(int i=0;i<3;i++) forceT[i] = forcesTan[index][i];
    return forceT;
}

// Capillary force
double* particleCloud::forceCap(int index) 
{
    double *forceC = new double[3];
    for(int i=0;i<3;i++) forceC[i] = forcesCap[index][i];
    return forceC;
}

// Viscous force
double* particleCloud::forceVis(int index) 
{
    double *forceV = new double[3];
    for(int i=0;i<3;i++) forceV[i] = forcesVis[index][i];
    return forceV;
}

double particleCloud::density(int index) 
{
  	double density;
  	density = densities[index];
  	return density; 
}

// Public member functions
void particleCloud::allocateVariables()
{
	positions = new double*[numP()];
	for(int i = 0; i < numP(); ++i) {
	    positions[i] = new double[3];
	}

	velocities = new double*[numP()];
	for(int i = 0; i < numP(); ++i) {
	    velocities[i] = new double[3];
	}
	
	omegas = new double*[numP()];
	for(int i = 0; i < numP(); ++i) {
	    omegas[i] = new double[3];
	}	
	
	liquidOns = new double[numP()];
	rads = new double[numP()];
		
	forces = new double*[numP()];
	for(int i = 0; i < numP(); ++i) 
	{
	    forces[i] = new double[3];
	}

	forcesTan = new double*[numP()];
	for(int i = 0; i < numP(); ++i) 
	{
	    forcesTan[i] = new double[3];
	}	

	forcesCap = new double*[numP()];
	for(int i = 0; i < numP(); ++i) 
	{
	    forcesCap[i] = new double[3];
	}
	
	forcesVis = new double*[numP()];
	for(int i = 0; i < numP(); ++i) 
	{
	    forcesVis[i] = new double[3];
	}	

	densities = new double[numP()];	
}

void particleCloud::setPos(double** pos)
{		
	for(int index = 0;index <  numP(); ++index)
    {
		for(int i=0;i<3;i++){
            positions[index][i] = pos[index][i];
        }
    }
	
}

void particleCloud::setVel(double** vel, double** ome)
{		
    for(int index = 0;index <  numP(); ++index)
    {
        for(int i=0;i<3;i++){
            velocities[index][i] = vel[index][i];
        }
    }

    for(int index = 0;index <  numP(); ++index)
    {
        for(int i=0;i<3;i++){
            omegas[index][i] = ome[index][i];
        }
    }	
	
}

void particleCloud::setLiquid(double* liq)
{
    for(int index = 0;index <  numP(); ++index)
    {
    	liquidOns[index] = liq[index];
    }
}

void particleCloud::setRadius(double* rad)
{
    for(int index = 0;index <  numP(); ++index)
    {
    	rads[index] = rad[index];
    }
}

void particleCloud::setForce(double** fcoll, double** ftan, double** fcap, double** fvis)
{			
    for(int index = 0;index <  numP(); ++index)
    {
        for(int i=0;i<3;i++){
            forces[index][i] = fcoll[index][i];
        }
    }

    for(int index = 0;index <  numP(); ++index)
    {
        for(int i=0;i<3;i++){
            forcesTan[index][i] = ftan[index][i];
        }
    }

    for(int index = 0;index <  numP(); ++index)
    {
        for(int i=0;i<3;i++){
            forcesCap[index][i] = fcap[index][i];
        }
    }
	
    for(int index = 0;index <  numP(); ++index)
    {
        for(int i=0;i<3;i++){
            forcesVis[index][i] = fvis[index][i];
        }
    }
				
}

void particleCloud::setDensity(double* dens)
{
    for(int index = 0;index <  numP(); ++index)
    {
    	densities[index] = dens[index];
    }
}


