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

	
      //evaluate the neural network model for the current point
      std::vector<float> Yresponse = YnnModel_.compute_output( &YdataChunk );
      Ytarget = Yresponse.at(0);
      if( Ytarget < Ymin ) Ytarget = Ymin;
      if( Ytarget > Ymax ) Ytarget = Ymax;   
      //print one point to log
      if(iCell==215)
	{
	  Info <<"iCell:"<<iCell<<endl;
	  Info << "Yfeature values:" << endl;
	  for(int i =0; i< YNFeatures; ++i) Info << Yfeatures[i] << " ";
	  Info << "Yresponse values:" << Yresponse.at(0) << endl;
	}
      if(iCell==11)
	{
	  Info <<"iCell:"<<iCell<<endl;
	  Info << "Yfeature values:" << endl;
	  for(int i =0; i< YNFeatures; ++i) Info << Yfeatures[i] << " ";
	  Info << "Yresponse values:" << Yresponse.at(0) << endl;
	}
      if(iCell==411)
	{
	  Info <<"iCell:"<<iCell<<endl;
	  Info << "Yfeature values:" << endl;
	  for(int i =0; i< YNFeatures; ++i) Info << Yfeatures[i] << " ";
	  Info << "Yresponse values:" << Yresponse.at(0) << endl;
	}	    
    }
  Ysgs[iCell].x() = 0.0;	
  Ysgs[iCell].y() = 0.0;
  if (solidVolumeFraction[iCell]<1e-3){
    Ysgs[iCell].z()=0.0;
  }
  else{
  Ysgs[iCell].z() = Ytarget;
  }	
}
