    //Create an implicit cylindrical porous zone not handled with fvOptions
    cylPorousZones pZones(mesh);
    if(pZones.size()>0)
    {
        Info << "Found cylindricalPorousZones" << endl;
        bool calcUCyl = false;
        if(pZones[0].dict().found("calcUCyl")) calcUCyl = true;

        vector axis(0, 0, 0);
        pZones[0].dict().readIfPresent("cylAxis",axis) ;
        if(mag(axis)<1) Info << "Warning: no valid axis supplied!" << endl;
        vector origin = pZones[0].origin() ;
        scalar omega = pZones[0].omega() ;
        Info << "pZones.origin() / axis() / omega() / calcUCyl " << origin << tab << axis  << tab << omega << tab << calcUCyl << endl;
    }
    else Info << "No cylindricalPorousZones" << endl;
    



