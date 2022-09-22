keras::DataChunkFlat DFXdataChunk( DFNFeatures, 0.0 );
keras::DataChunkFlat DFYdataChunk( DFNFeatures, 0.0 );
keras::DataChunkFlat DFZdataChunk( DFNFeatures, 0.0 );
std::vector<float>& DFXfeatures = DFXdataChunk.get_1d_rw();
std::vector<float>& DFYfeatures = DFYdataChunk.get_1d_rw();
std::vector<float>& DFZfeatures = DFZdataChunk.get_1d_rw();

forAll( driftFlux_, iCell )
{

  float DFXtarget = 1.0e-6;
  float DFYtarget = 1.0e-6;
  float DFZtarget = 1.0e-6;
  unsigned int index = 0;
  DFXfeatures.clear();
  DFYfeatures.clear();
  DFZfeatures.clear();
   //solid volume fraction at the cell 
   DFXfeatures.push_back( ( solidVolumeFraction[iCell]/alphaMax - DFXmeans_[index])/DFXvars_[index] );
   DFYfeatures.push_back( ( solidVolumeFraction[iCell]/alphaMax - DFYmeans_[index])/DFYvars_[index] );
   DFZfeatures.push_back( ( solidVolumeFraction[iCell]/alphaMax - DFZmeans_[index])/DFZvars_[index] );
   ++index;

   DFXfeatures.push_back((gradP[iCell].x()/9.8/rhoParticle_-DFXmeans_[index])/DFXvars_[index]);
   DFYfeatures.push_back((gradP[iCell].y()/9.8/rhoParticle_-DFYmeans_[index])/DFYvars_[index]);
   DFZfeatures.push_back((gradP[iCell].z()/9.8/rhoParticle_-DFZmeans_[index])/DFZvars_[index]);
   ++index;

   //uSlip.z at the center 
   DFXfeatures.push_back((uSlip[iCell].x()/uTerminal-DFXmeans_[index])/DFXvars_[index]);
   DFYfeatures.push_back((uSlip[iCell].y()/uTerminal-DFYmeans_[index])/DFYvars_[index]);
   DFZfeatures.push_back((uSlip[iCell].z()/uTerminal-DFZmeans_[index])/DFZvars_[index]);
   ++index;
   //evaluate the neural network model for the current point
   std::vector<float> DFXresponse = DFXnnModel_.compute_output( &DFXdataChunk );
   std::vector<float> DFYresponse = DFYnnModel_.compute_output( &DFYdataChunk );
   std::vector<float> DFZresponse = DFZnnModel_.compute_output( &DFZdataChunk );
   DFXtarget = DFXresponse.at(0);
   DFYtarget = DFYresponse.at(0);
   DFZtarget = DFZresponse.at(0);
   /*
   if(iCell==215)
     {
	Info<<"icell:"<<iCell<<endl;
	Info << "DFfeature values:" << endl;
	for(int i =0; i<DFNFeatures; ++i) Info << DFfeatures[i] << " ";
	Info << "DFresponse values:" << DFresponse.at(0) << endl;
     }
   */
   driftFlux_[iCell].x() = DFXtarget*0.64*0.219;
   driftFlux_[iCell].y() = DFYtarget*0.64*0.219;
   driftFlux_[iCell].z() = DFZtarget*0.64*0.219;

}
