{

//features
//uSlip = pair_.continuous().U() - pair_.dispersed().U();
//uSlipVdrift = pair_.continuous().U() - pair_.dispersed().U()+driftFlux_;
//gradAlphaP = fvc::grad( solidVolumeFraction );

stencil.collectData( solidVolumeFraction, alphaList_ );
stencil.collectData( driftFlux_, driftFluxList_ );
stencil.collectData( gradAlphaP, gradAlphaPList_ );
stencil.collectData( uSlip, uSlipList_ );

Info << "In Ysgs_z.H" << endl;    
    
keras::DataChunkFlat YdataChunk( YNFeatures, 0.0 );
    
std::vector<float>& Yfeatures = YdataChunk.get_1d_rw();
    
forAll( alphaList_, iCell )
{
        
  List<scalar>& alphaSubList_ = alphaList_[iCell];
  List<vector>& driftFluxSubList_ = driftFluxList_[iCell];
  List<vector>& uSlipSubList_ = uSlipList_[iCell];
  List<vector>& gradAlphaPSubList_ = gradAlphaPList_[iCell];
  //List<vector>& uSlipVdriftSubList_ = uSlipVdriftList_[iCell];

  //collect features from the stencil
  float Ytarget = 1.0e-6;
  float gradYtarget = 0.0;
  float driftFluxEq_ = 0.0;
  	
  if( alphaSubList_.size() >= 27 )
  {

    /* INSERT features */

    unsigned int index = 0;

    Yfeatures.clear();


    //solid volume fraction at the cell 
    Yfeatures.push_back( ( alphaSubList_[0]/alphaMax - Ymeans_[index])/Yvars_[index] );
    ++index;

    //gradAlphaP at center cell
    Yfeatures.push_back((gradAlphaPSubList_[0].x()*lStar-Ymeans_[index])/Yvars_[index]);
    ++index;
    Yfeatures.push_back((gradAlphaPSubList_[0].y()*lStar-Ymeans_[index])/Yvars_[index]);
    ++index;
    Yfeatures.push_back((gradAlphaPSubList_[0].z()*lStar-Ymeans_[index])/Yvars_[index]);
    ++index;	    


    //Vdrift.z at the center 
    Yfeatures.push_back((driftFluxSubList_[0].z()/uTerminal/alphaMax-Ymeans_[index])/Yvars_[index]);
    ++index;

    //uSlip.z at the center 
    Yfeatures.push_back((uSlipSubList_[0].z()/uTerminal-Ymeans_[index])/Yvars_[index]);
    ++index;    
    
    //solve X + Y == 0
    
    driftFluxEq_ = 0.0;
    scalar driftFluxDelta( 10 );
    int iteration = 0;
    int maxIter = 10;
    
    while( driftFluxDelta >= 1e-7 && iteration < maxIter )
    {
       
	++iteration;
	
	float driftFluxOld = driftFluxEq_;
	
	Yfeatures[4] = (driftFluxEq_/uTerminal/alphaMax-Ymeans_[4])/Yvars_[4]; 
	
    	std::vector<float> Yresponse = YnnModel_.compute_output( &YdataChunk );
        Ytarget = Yresponse.at(0);	
	
	float ssXY = Xsgs[iCell].z()/(1.0-solidVolumeFraction[iCell]) + Ytarget;
	
	//Info<<ssXY<<"    ";
	
	#include "gradYsgs.H"
	
	if( gradYtarget < 1e-7 && gradYtarget > -1e-7 )
	{
	   driftFluxEq_ = 0.0;
	   Info<< "Warning: Newton Iteration Diverged! ( "<< gradYtarget<<" "<<ssXY<<" )"<<endl;
	   break;
	}
	
	driftFluxEq_ = driftFluxOld - 5e-2 * ssXY/gradYtarget;
	
	while( driftFluxEq_ * driftFluxOld < 0 )
	{
	   driftFluxEq_ = driftFluxOld * 0.5 + driftFluxEq_ * 0.5;
	}
	
	driftFluxDelta = mag( scalar( driftFluxEq_ - driftFluxOld ) );
	
    }
    
    
    //Info<< "DriftFluxEq = " << driftFluxEq_ << " ( " << iteration << " )" << endl;
    
        
  }
  
  
  if (solidVolumeFraction[iCell]<1e-3)
  {
    driftFlux_[iCell].x() = 0.0;
    driftFlux_[iCell].y() = 0.0;
    driftFlux_[iCell].z() = 0.0;    
  }
  else
  {
    driftFlux_[iCell].x() = 0.0;
    driftFlux_[iCell].y() = 0.0;
    driftFlux_[iCell].z() = driftFluxEq_;
  }	
  
}

}
