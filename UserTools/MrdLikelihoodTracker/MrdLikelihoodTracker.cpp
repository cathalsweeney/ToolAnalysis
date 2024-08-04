#include "MrdLikelihoodTracker.h"


// Return the (x,y) coordinates of a line at a specified z-value,
// given the parameters of the line
std::pair<double, double> CalcCoordsAtZ(double z_val,   double x_start, double y_start,
                                    double z_start, double theta,   double phi)
{

  double sin_phi = std::sin(phi*(M_PI/180.) );
  double cos_phi = std::cos(phi*(M_PI/180.) );
  double tan_theta = std::tan(theta*(M_PI/180.) );
  
  double x_val = x_start + (z_val - z_start)*tan_theta*cos_phi;
  double y_val = y_start + (z_val - z_start)*tan_theta*sin_phi;
//  std::cout << "{ " << x_val << ", " << y_val << ", " << z_val << "} \n";

  return {x_val, y_val};
}

MrdLikelihoodTracker::MrdLikelihoodTracker():Tool(){}


bool MrdLikelihoodTracker::Initialise(std::string configfile, DataModel &data)
{

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  if(!m_variables.Get("Verbose",fVerbose)){
    fVerbose = 0;
    Log("MrdLikelihoodTracker tool: Did not find \"Verbose\" in config file. Using default value of 0",v_message,fVerbose); 
  }
  else if(fVerbose < 0){
    fVerbose = 0;
    Log("MrdLikelihoodTracker tool: \"Verbose\" cannot be < 0. Using default value of 0",v_error,fVerbose);
  }

  
  if(!m_variables.Get("Event",fEventNumber)){
    Log("MrdLikelihoodTracker tool: Did not find \"Event\" in config file. Using default value of 1",v_message,fVerbose); 
    fEventNumber = 1;
  }
  else{
    std::cout << "Hi " << fEventNumber << "\n";
    if(fEventNumber < 0){
      Log("MrdLikelihoodTracker tool: \"Event\" cannot be < 1. Quitting",v_error,fVerbose);
      return false;
    }
  }


  if(!m_variables.Get("StartX",fStartX)){
    Log("MrdLikelihoodTracker tool: Did not find \"StartX\" in config file. Quitting",v_message,fVerbose); 
    return false;
  }
  if(!m_variables.Get("StartY",fStartY)){
    Log("MrdLikelihoodTracker tool: Did not find \"StartY\" in config file. Quitting",v_message,fVerbose); 
    return false;
  }
//  if(!m_variables.Get("StartZ",fStartZ)){
//    Log("MrdLikelihoodTracker tool: Did not find \"StartZ\" in config file. Quitting",v_message,fVerbose); 
//    return false;
//  }

  
  if(!m_variables.Get("Theta",fTheta)){
    Log("MrdLikelihoodTracker tool: Did not find \"Theta\" in config file. Quitting",v_message,fVerbose); 
    return false;
  }  
  if(!m_variables.Get("Phi",fPhi)){
    Log("MrdLikelihoodTracker tool: Did not find \"Phi\" in config file. Quitting",v_message,fVerbose); 
    return false;
  }

//  hist = new TH1D("hist", "hist", 200, 3., 5.);
  h2 = new TH2D("hist", "hist", nBinsX,0.,60., nBinsY,0.,359.);

  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",fGeom);
  int n_mrd_pmts = fGeom->GetNumDetectorsInSet("MRD");

  std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = fGeom->GetDetectors();
  for(std::map<unsigned long,Detector*>::iterator it  = Detectors->at("MRD").begin();
      it != Detectors->at("MRD").end();
      ++it){
    unsigned long detkey = it->first; // I'm pretty sure detkey == chankey
    fAllMrdChankeys.push_back(detkey);
  }
//  delete Detectors;
  
  return true;
}


bool MrdLikelihoodTracker::Execute()
{
  Log("MrdLikelihoodTracker tool: Executing.",v_message,fVerbose);

  if(!m_data->Stores.count("ANNIEEvent") ){ 
    Log("MrdLikelihoodTracker tool: No ANNIEEvent store!",v_error,fVerbose); 
    return false;
  }

  int this_event;
  int get_ok = 0;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("MCEventNum",this_event);
  if(not get_ok){ Log("MrdLikelihoodTracker tool: Error retrieving MCEventNum, true from ANNIEEvent!",v_error,fVerbose); return false; }

  if(this_event != 57) return true; 
  std::cout << "this_event is " << this_event << "\n";
  
  std::vector<std::vector<int>> MrdTimeClusters;
  get_ok = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
  if (not get_ok) { Log("MrdLikelihoodTracker Tool: Error retrieving MrdTimeClusters map from CStore, did you run TimeClustering beforehand?",v_error,fVerbose); return false; }

  int nClust = MrdTimeClusters.size();
  std::cout << "Number of MrdTimeClusters " << nClust << "\n";
  if(nClust != 1) return false; // TODO maybe add a log message
  std::cout << "Size of MrdTimeClusters[0] " << MrdTimeClusters[0].size() << "\n";


  //get_ok = m_data->CStore.Get("MrdDigitTimes",MrdDigitTimes);
//  if (not get_ok) { Log("EventDisplay Tool: Error retrieving MrdDigitTimes map from CStore, did you run TimeClustering beforehand?",v_error,verbose); return false; }
  get_ok = m_data->CStore.Get("MrdDigitChankeys",mrddigitchankeysthisevent);
  if (not get_ok) { Log("EventDisplay Tool: Error retrieving MrdDigitChankeys, did you run TimeClustering beforehand",v_error,fVerbose); return false;}

  std::vector<int> hitmrd_chankeys;
  
  for(int digit_value : MrdTimeClusters[0]){
//    std::cout << digit_value << "\n";
    unsigned long chankey = mrddigitchankeysthisevent.at(digit_value);
//    Detector *thedetector = fGeom->ChannelToDetector(chankey);
//    unsigned long detkey = thedetector->GetDetectorID();
//    if (thedetector->GetDetectorElement()!="MRD") continue;
//    double mrdtimes=MrdDigitTimes.at(digit_value);
    hitmrd_chankeys.push_back(chankey);    
//    std::cout << chankey << "\n";
  }


  for(int iBinX=1; iBinX<=nBinsX; iBinX++){
    double binX = h2->GetXaxis()->GetBinCenter(iBinX);
    fTheta = binX;
    for(int iBinY=1; iBinY<=nBinsY; iBinY++){
      double binY = h2->GetYaxis()->GetBinCenter(iBinY);
      fPhi = binY;

      //~~~~~~~~~~~~~~~~~~~~~~~~~
      FillCoordsAtZ();
      FillPaddleProbs();
      
      double prob = 1.;
//      std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = fGeom->GetDetectors();
//      for(std::map<unsigned long,Detector*>::iterator it  = Detectors->at("MRD").begin();
//          it != Detectors->at("MRD").end();
//          ++it){
//        unsigned long detkey = it->first; // I'm pretty sure detkey == chankey
//        int chankey = detkey; // redundant, but makes code more readable
      for(int chankey : fAllMrdChankeys){

        bool isHit = false;
        double val = 1.;
        if(std::find(hitmrd_chankeys.begin(), hitmrd_chankeys.end(), chankey) != hitmrd_chankeys.end()){
          val = fPaddleProbs.at(chankey); 
          isHit = true;
        }
        else{
          val = ( 1.- fPaddleProbs.at(chankey) );
        }    
          prob *= val;

//          std::cout << val << "\n";
//          if(val == 1. || val < 0.){
//            std::cout << isHit << " : " << fPaddleProbs.at(chankey) << "\n";
//          }

      }// end loop over all MRD channels
//      std::cout << "prob is " << prob << "\n";
      //~~~~~~~~~~~~~~~~~~~~~~~~~
      h2->SetBinContent(iBinX, iBinY, prob);
    }//end loop over y-bins
  }//end loop over x-bins
  
  return true;
}

void MrdLikelihoodTracker::FillCoordsAtZ()
{
  fCoordsAtZ.clear(); // TODO might be unecessary?
  for(double z_mid : fZ_midpoints){
    fCoordsAtZ[z_mid] = CalcCoordsAtZ(z_mid, fStartX, fStartY, fZ_start, fTheta, fPhi);
  }
  
}
  
void MrdLikelihoodTracker::FillPaddleProbs()
{
  fPaddleProbs.clear(); // TODO might be unecessary?


//  std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = fGeom->GetDetectors();
//  for(std::map<unsigned long,Detector*>::iterator it  = Detectors->at("MRD").begin();
//                                                    it != Detectors->at("MRD").end();
//                                                  ++it){
//  Detector* amrdpmt = it->second;
//    unsigned long detkey = it->first;
//    unsigned long chankey = amrdpmt->GetChannels()->begin()->first;
  
  for(int chankey : fAllMrdChankeys){

    Paddle *mrdpaddle = (Paddle*) fGeom->GetDetectorPaddle(chankey);

//    if(chankey != detkey) std::cout << "HELLO ( " << chankey << ", " << detkey << ") \n";
    
    double x_min = mrdpaddle->GetXmin();
    double x_max = mrdpaddle->GetXmax();
    double x_mid = (x_min + x_max)/ 2.0; // TODO can i just use paddle->origin instead?
    double x_extent = std::abs(x_max - x_min); // TODO theres prob already a member for this
         
    double y_min = mrdpaddle->GetYmin();
    double y_max = mrdpaddle->GetYmax();
    double y_mid = (y_min + y_max)/ 2.0; // TODO can i just use paddle->origin instead?
    double y_extent = std::abs(y_max - y_min); // TODO theres prob already a member for this
    
    double z_min = mrdpaddle->GetZmin();
    double z_max = mrdpaddle->GetZmax();
    double z_mid = (z_min + z_max)/2.0;

    std::pair<double, double> coords = fCoordsAtZ[z_mid];

    //TODO bad name
    double chi2_x = (coords.first - x_mid)*(coords.first - x_mid) / (x_extent * x_extent);
    double chi2_y = (coords.second - y_mid)*(coords.second - y_mid) / (y_extent * y_extent);


    double prob = ( 1./(2*M_PI*x_extent*y_extent) ) * exp(-chi2_x/2. - chi2_y/2.);
    if(prob > 1.){
      std::cout << "~~~~ \n";
      std::cout << prob << "\n";
      std::cout << chi2_x << " : (" << coords.first << " - " << x_mid << ") / " << x_extent << "\n";
      std::cout << chi2_y << " : (" << coords.second << " - " << y_mid << ") / " << y_extent << "\n"; 
    }
    
    fPaddleProbs[chankey] = prob;
//    std::cout << chankey << ", " << prob << "\n";
    
    
//    int orientation = mrdpaddle->GetOrientation();    //0 is horizontal, 1 is vertical
//    int half = mrdpaddle->GetHalf();                  //0 or 1
//    int side = mrdpaddle->GetSide();


//    std::cout << "( " << z_min << ", " << z_max << ") \n";
//    std::cout << z_mid << "\n";
//    hist->Fill(z_mid);
//    if(std::find(fZ_midpoints.begin(), fZ_midpoints.end(), z_mid) == fZ_midpoints.end()){
//      fZ_midpoints.push_back(z_mid);
//    }

  }// end loop over MRD paddles
  

  
}


bool MrdLikelihoodTracker::Finalise()
{

//  std::cout << "Midpoints \n";
//  for(double m : fZ_midpoints) std::cout << m << "\n";

  TCanvas* c = new TCanvas();
  c->SetLogz(true);
  c->SetRightMargin(0.13);
  
  h2->SetStats(0);
  h2->GetZaxis()->SetRangeUser(1e-25, 1e-22);
  h2->Draw("colz");
  c->SaveAs("bar.png");

  delete h2;
  h2 = nullptr;

  delete c;
  c = nullptr;
  
  return true;
}
