
keras::DataChunkFlat XdataChunk( XNFeatures, 0.0 );
    
std::vector<float>& Xfeatures = XdataChunk.get_1d_rw();
    
forAll( solidVolumeFraction(), iCell )
{
  //collect features from the stencil
  float Xtarget = 1.0e-6;
	    
      /* INSERT features */
      
  unsigned int index = 0;
  Xfeatures.clear();
	    
  //solid volume fraction at the cell 
  Xfeatures.push_back( ( solidVolumeFraction()[iCell]/alphaMax - Xmeans_[index])/Xvars_[index] );
  ++index;

  //gradAlphaP at the center
  Xfeatures.push_back((gradAlphaP()[iCell].x()*lStar-Xmeans_[index])/Xvars_[index]);
  ++index;
  Xfeatures.push_back((gradAlphaP()[iCell].y()*lStar-Xmeans_[index])/Xvars_[index]);
  ++index;
  Xfeatures.push_back((gradAlphaP()[iCell].z()*lStar-Xmeans_[index])/Xvars_[index]);
  ++index;

  //gradP at the center 
  Xfeatures.push_back((gradP()[iCell].z()/9.8/rhoParticle_-Xmeans_[index])/Xvars_[index]);
  ++index;


  //gradGradPz at the center
  Xfeatures.push_back((gradGradPz()[iCell].x()/9.8/rhoParticle_*lStar-Xmeans_[index])/Xvars_[index]);
  ++index;
  Xfeatures.push_back((gradGradPz()[iCell].y()/9.8/rhoParticle_*lStar-Xmeans_[index])/Xvars_[index]);
  ++index;
  Xfeatures.push_back((gradGradPz()[iCell].z()/9.8/rhoParticle_*lStar-Xmeans_[index])/Xvars_[index]);
  ++index;

  //laplaceGradPgZ at the center
  Xfeatures.push_back((laplaceGradPgZ()[iCell]/9.8/rhoParticle_*lStar*lStar-Xmeans_[index])/Xvars_[index]);
      
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
  /*
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
  */
  Xsgs[iCell].x() = 0.0;	
  Xsgs[iCell].y() = 0.0;
  if(solidVolumeFraction()[iCell]<1e-3){
    Xsgs[iCell].z()=0.0;
  }else{
    Xsgs[iCell].z() = Xtarget;
  }
	
}
