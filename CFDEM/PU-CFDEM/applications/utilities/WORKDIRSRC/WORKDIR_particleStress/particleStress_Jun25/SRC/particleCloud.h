class particleCloud 
{
    int width, height;
	double partR;
	int partID;
	int nP;
    mutable double **positions;
    mutable double **velocities;
	mutable double **omegas;
	mutable	double *liquidOns;
	mutable	double *rads;
	mutable double **forces;
	mutable double **forcesTan;
	mutable double **forcesCap;	
	mutable double **forcesVis;		
	mutable double *densities;

  public:
    void set_values (int,int);

	int setNumberOfParticles (int);
	
    double* position(int);
    double* velocity(int);
    double* omega(int);
	double liquidOn(int);
	double radius(int);
	double* force(int);
	double* forceTan(int);
	double* forceCap(int);
	double* forceVis(int);	
	double density(int);
		
    int area() {return width*height;}
	int numP() {return nP;}
	
    void setPos(double **);
	void setVel(double **, double **);	
	void setLiquid(double *);
	void setRadius(double *);
	void setForce(double **, double **, double **, double **);	
	void setDensity(double *);
	
	void allocateVariables();	
		
};