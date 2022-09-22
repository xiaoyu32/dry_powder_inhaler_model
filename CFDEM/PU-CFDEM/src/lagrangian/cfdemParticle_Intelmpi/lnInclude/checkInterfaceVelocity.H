//check the interface velocity of particle
    int count=0;
    scalar summag=0; 
    vector Usur(0,0,0);
    vector Usum(0,0,0);
    scalar ubar=0;
    int numprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    OFstream* interfaceProps_;
    interfaceProps_ =  new OFstream("interfaceProps.dat");
    *interfaceProps_  << "#Time interface no-slip checking" << endl;

    for(int index=0; index< numberOfParticles(); index++)
    {
            for(int subCell=0;subCell<voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                cellI = cellIDs()[index][subCell];
                Info << "subcell=" << cellIDs()[index][subCell] << endl;
                if (cellI >= 0 &&voidfractions_[index][subCell]>0&&voidfractions_[index][subCell]<1)
                {
                     count++; 
                     //Usur = voidfractions_[index][subCell]*U[cellI] - voidfractions_[index][subCell]*uParticle;
                     Usur = voidfractions_[index][subCell]*U[cellI] - uParticle;
                     summag += mag(Usur);
                     //Info << "erf=" << erf(mag(Usur))<< endl;
                }
            }
    }

    ubar = summag/count;
    printf("process %d \n",me);
    Info << "ubar=" << ubar << endl;
    Info << "count=" << count << endl;

    if(Pstream::master()) //Write only if master
    {
        *interfaceProps_ << setw(IOstream::defaultPrecision() + 6) 
                << U.mesh().time().value() << "   " 
                << ubar    << "   "
                << endl;
    }

