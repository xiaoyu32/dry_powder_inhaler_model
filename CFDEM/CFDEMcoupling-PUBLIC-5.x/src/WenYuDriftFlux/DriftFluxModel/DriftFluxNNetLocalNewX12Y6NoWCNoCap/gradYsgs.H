{
   
   
   
   float deltaVd = 1e-2;
   float vdOrig = Yfeatures[4];
   
   Yfeatures[4] = vdOrig + deltaVd;
   Yresponse = YnnModel_.compute_output( &YdataChunk );   
   float Ytarget2 =  Yresponse.at(0);

   Yfeatures[4] = vdOrig - deltaVd;
   Yresponse = YnnModel_.compute_output( &YdataChunk );   
   float Ytarget1 =  Yresponse.at(0);
   
   gradYtarget = ( Ytarget2 - Ytarget1 )/( 2.0 * deltaVd ); 
   
   //Info<< "Y gradient = "<<gradYtarget<<"  "<<Ytarget2<<"  "<<Ytarget1<< endl;
   
}
