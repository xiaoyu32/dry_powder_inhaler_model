//const extendedCentredCellToCellStencil& stencil = this->stencil();

//gradP = fvc::grad( pair_.continuous().thermo().p());
//gradAlphaP = fvc::grad( solidVolumeFraction );


//features
stencil.collectData( solidVolumeFraction, alphaList_ );
stencil.collectData( gradP, gradPList_ );
stencil.collectData( gradAlphaP, gradAlphaPList_ );
    
//stencil.collectData( mesh_.C(), uSolidList_ );
    
keras::DataChunkFlat XdataChunk( XNFeatures, 0.0 );
keras::DataChunkFlat XBdataChunk( XBNFeatures,0.0 );
    
std::vector<float>& Xfeatures = XdataChunk.get_1d_rw();
std::vector<float>& XBfeatures = XBdataChunk.get_1d_rw();    
forAll( alphaList_, iCell )
{
        
  List<scalar>& alphaSubList_ = alphaList_[iCell];
  List<vector>& gradPSubList_ = gradPList_[iCell];
  List<vector>& gradAlphaPSubList_ = gradAlphaPList_[iCell];
	
  //collect features from the stencil
  float Xtarget = 1.0e-6;
	
  if( alphaSubList_.size() >= 27 )
    {
	    
      /* INSERT features */
      
      unsigned int index = 0;
      Xfeatures.clear();
      
      //solid volume fraction at z-1 
      Xfeatures.push_back( ( alphaSubList_[11]/alphaMax - Xmeans_[index])/Xvars_[index] );
      ++index;
	    
      //gradP at z-1 
      Xfeatures.push_back((gradPSubList_[11].x()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradPSubList_[11].y()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradPSubList_[11].z()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
      ++index;

      //gradAlphaP at z-1
      Xfeatures.push_back((gradAlphaPSubList_[11].x()*lStar-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradAlphaPSubList_[11].y()*lStar-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradAlphaPSubList_[11].z()*lStar-Xmeans_[index])/Xvars_[index]);
      ++index;
      
	    
      //solid volume fraction at the cell 
      Xfeatures.push_back( ( alphaSubList_[0]/alphaMax - Xmeans_[index])/Xvars_[index] );
      ++index;
	    
      //gradP at the center 
      Xfeatures.push_back((gradPSubList_[0].x()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradPSubList_[0].y()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradPSubList_[0].z()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
      ++index;

      //gradAlphaP at the center
      Xfeatures.push_back((gradAlphaPSubList_[0].x()*lStar-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradAlphaPSubList_[0].y()*lStar-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradAlphaPSubList_[0].z()*lStar-Xmeans_[index])/Xvars_[index]);
      ++index;
      	    
      //solid volume fraction at z+1 
      Xfeatures.push_back( ( alphaSubList_[22]/alphaMax - Xmeans_[index])/Xvars_[index] );
      ++index;
	    
      //gradP at z+1 
      Xfeatures.push_back((gradPSubList_[22].x()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradPSubList_[22].y()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradPSubList_[22].z()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
      ++index;

      //gradAlphaP at z+1
      Xfeatures.push_back((gradAlphaPSubList_[22].x()*lStar-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradAlphaPSubList_[22].y()*lStar-Xmeans_[index])/Xvars_[index]);
      ++index;
      Xfeatures.push_back((gradAlphaPSubList_[22].z()*lStar-Xmeans_[index])/Xvars_[index]);
      ++index;
      
      //evaluate the neural network model for the current point
      std::vector<float> Xresponse = XnnModel_.compute_output( &XdataChunk );
      Xtarget = Xresponse.at(0);
      if( Xtarget < Xmin ) Xtarget=Xmin;
      if( Xtarget > Xmax ) Xtarget=Xmax;   
      //print one point to log
      if(iCell==215)
	{
	  Info<<"icell:"<<iCell<<endl;
	  Info << "Xfeature values:" << endl;
	  for(int i =0; i<XNFeatures; ++i) Info << Xfeatures[i] << " ";
	  Info << "Xresponse values:" << Xresponse.at(0) << endl;
	}
      if(iCell==11)
	{
	  Info<<"icell:"<<iCell<<endl;
	  Info << "Xfeature values:" << endl;
	  for(int i =0; i<XNFeatures; ++i) Info << Xfeatures[i] << " ";
	  Info << "Xresponse values:" << Xresponse.at(0) << endl;
	}
      if(iCell==411)
	{
	  Info<<"icell:"<<iCell<<endl;
	  Info << "Xfeature values:" << endl;
	  for(int i =0; i<XNFeatures; ++i) Info << Xfeatures[i] << " ";
	  Info << "Xresponse values:" << Xresponse.at(0) << endl;
	}    
    }
  Xsgs[iCell].x() = 0.0;	
  Xsgs[iCell].y() = 0.0;
  if(solidVolumeFraction[iCell]<1e-3){
    Xsgs[iCell].z()=0.0;
  }else{
    Xsgs[iCell].z() = Xtarget;
  }	
}
