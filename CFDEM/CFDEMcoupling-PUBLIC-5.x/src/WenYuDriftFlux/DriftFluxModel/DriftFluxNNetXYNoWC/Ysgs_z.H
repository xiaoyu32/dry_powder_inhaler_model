//features
//uSlip = pair_.continuous().U() - pair_.dispersed().U();
//uSlipVdrift = pair_.continuous().U() - pair_.dispersed().U()+driftFlux_;
//gradAlphaP = fvc::grad( solidVolumeFraction );
stencil.collectData( solidVolumeFraction, alphaList_ );
stencil.collectData( gradP, gradPList_ );
stencil.collectData( gradAlphaP, gradAlphaPList_ );


Info << "In Ysgs_z.H" << endl;    
    
keras::DataChunkFlat YdataChunk( YNFeatures, 0.0 );
    
std::vector<float>& Yfeatures = YdataChunk.get_1d_rw();
    
forAll( alphaList_, iCell )
{
    List<scalar>& alphaSubList_ = alphaList_[iCell];
  List<vector>& gradPSubList_ = gradPList_[iCell];
  List<vector>& gradAlphaPSubList_ = gradAlphaPList_[iCell];      


  //collect features from the stencil
  float Ytarget = 1.0e-6;
	
  if( alphaSubList_.size() >=0)//>= 27 )
    {

      /* INSERT features */
      
      unsigned int index = 0;
      Yfeatures.clear();
	    
      //solid volume fraction at the cell 
      Yfeatures.push_back( ( alphaSubList_[0]/alphaMax - Ymeans_[index])/Yvars_[index] );
      ++index;

      //gradAlphaP at the center
      Yfeatures.push_back((gradAlphaPSubList_[0].x()*lStar-Ymeans_[index])/Yvars_[index]);
      ++index;
      Yfeatures.push_back((gradAlphaPSubList_[0].y()*lStar-Ymeans_[index])/Yvars_[index]);
      ++index;
      Yfeatures.push_back((gradAlphaPSubList_[0].z()*lStar-Ymeans_[index])/Yvars_[index]);
      ++index;

	    
      //gradP at the center 
      Yfeatures.push_back((gradPSubList_[0].x()/9.8/rhoParticle_-Ymeans_[index])/Yvars_[index]);
      ++index;
      Yfeatures.push_back((gradPSubList_[0].y()/9.8/rhoParticle_-Ymeans_[index])/Yvars_[index]);
      ++index;
      Yfeatures.push_back((gradPSubList_[0].z()/9.8/rhoParticle_-Ymeans_[index])/Yvars_[index]);
      ++index;



      //gradGradPz at the center
      Yfeatures.push_back((gradGradPz()[iCell].x()/9.8/rhoParticle_*lStar-Ymeans_[index])/Yvars_[index]);
      ++index;
      Yfeatures.push_back((gradGradPz()[iCell].y()/9.8/rhoParticle_*lStar-Ymeans_[index])/Yvars_[index]);
      ++index;
      Yfeatures.push_back((gradGradPz()[iCell].z()/9.8/rhoParticle_*lStar-Ymeans_[index])/Yvars_[index]);
      ++index;

      //laplaceGradPz at the center
      Yfeatures.push_back((laplaceGradPz()[iCell]/9.8/rhoParticle_*lStar*lStar-Ymeans_[index])/Yvars_[index]);
      ++index;
      Yfeatures.push_back((laplaceAlphaP()[iCell]*lStar*lStar-Ymeans_[index])/Yvars_[index]);
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
