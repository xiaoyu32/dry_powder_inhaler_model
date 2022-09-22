stencil.collectData( solidVolumeFraction, alphaList_ );
stencil.collectData( uSlip, uSlipList_ );

keras::DataChunkFlat DFdataChunk( DFNFeatures, 0.0 );
std::vector<float>& DFfeatures = DFdataChunk.get_1d_rw();

forAll( alphaList_, iCell )
{
  List<scalar>& alphaSubList_ = alphaList_[iCell];
  List<vector>& uSlipSubList_ = uSlipList_[iCell];

  //collect features from the stencil
  float DFtarget = 1.0e-6;

  if( alphaSubList_.size() >= 27 )
    {
	    
      /* INSERT features */
		    
      unsigned int index = 0;
	
      DFfeatures.clear();

      //(1)
      DFfeatures.push_back( ( alphaSubList_[13]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[7]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[4]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[21]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[16]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[10]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[1]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[24]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[18]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[14]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[8]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[5]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[22]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[0]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[11]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[2]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[25]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[19]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[9]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[15]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[6]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[23]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[17]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[12]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[3]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[26]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( alphaSubList_[20]/alphaMax - DFmeans_[index])/DFvars_[index] );
      ++index; 



      DFfeatures.push_back( ( uSlipSubList_[13].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[7].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[4].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[21].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[16].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[10].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[1].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[24].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[18].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[14].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[8].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[5].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[22].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[0].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[11].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[2].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[25].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[19].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[9].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[15].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[6].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[23].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[17].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[12].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[3].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[26].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;
      DFfeatures.push_back( ( uSlipSubList_[20].z()/uTerminal - DFmeans_[index])/DFvars_[index] );
      ++index;     

   //evaluate the neural network model for the current point
      std::vector<float> DFresponse = DFnnModel_.compute_output( &DFdataChunk );
      DFtarget = DFresponse.at(0);
      if(iCell==215)
	{
	  Info<<"icell:"<<iCell<<endl;
	  Info << "DFfeature values:" << endl;
	  for(int i =0; i<DFNFeatures; ++i) Info << DFfeatures[i] << " ";
	  Info << "DFresponse values:" << DFresponse.at(0) << endl;
	}
    }
   driftFlux_[iCell].x()=0.0;
   driftFlux_[iCell].y()=0.0;
   driftFlux_[iCell].z()=DFtarget*0.64*0.219;

}
