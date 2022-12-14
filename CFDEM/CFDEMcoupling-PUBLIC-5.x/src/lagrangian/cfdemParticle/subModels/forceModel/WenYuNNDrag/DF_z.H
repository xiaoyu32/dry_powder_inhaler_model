keras::DataChunkFlat DFdataChunk( DFNFeatures, 0.0 );
std::vector<float>& DFfeatures = DFdataChunk.get_1d_rw();

float DFtarget = 1.0e-6;
unsigned int index = 0;
DFfeatures.clear();

DFfeatures.push_back(( Resys_ -DFmeans_[index])/DFvars_[index]);
++index;

DFfeatures.push_back( (fs_-DFmeans_[index])/DFvars_[index]); 
++index;

DFfeatures.push_back( ( (1.0-voidfraction)/alphaMax - DFmeans_[index])/DFvars_[index] );
++index;

DFfeatures.push_back((gradP.z()/9.8/rho_s-DFmeans_[index])/DFvars_[index]);
++index;

DFfeatures.push_back((Ur.z()/uTerminal-DFmeans_[index])/DFvars_[index]);
++index;

std::vector<float> DFresponse = DFnnModel_.compute_output( &DFdataChunk );
DFtarget = DFresponse.at(0);

driftFlux = DFtarget*alphaMax*uTerminal;
