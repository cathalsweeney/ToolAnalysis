#include "MrdLikelihoodTracker.h"
#include "TMarker.h"


//..................................................................................
// Return the (x,y) coordinates of a line at a specified z-value,
// given the parameters of the line. 
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

//..................................................................................

MrdLikelihoodTracker::MrdLikelihoodTracker():Tool(){}
//..................................................................................

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

//  hist = new TH1D("hist", "hist", 200, 3., 5.);
//  h2 = new TH2D("hist", "hist", nBinsX,0.,60., nBinsY,0.,359.);
  h2_angle = new TH2D("", "", nBinsX,0.,90., nBinsY,0.,360.);
  h2_space = new TH2D("", "", nBinsX,-2.0,2.0, nBinsY,-2.0,2.0);
    

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

//..................................................................................

bool MrdLikelihoodTracker::Execute()
{  
  Log("MrdLikelihoodTracker tool: Executing.",v_message,fVerbose);

  // TODO could move these lines to a Reset() method
  fMinimizer = NULL;
  fFitStatus = 0;
  fFitVals.clear();
//  fStartX = fStartX_default;
//  fStartY = fStartY_default;
//  fTheta = fTheta_default;
//  fPhi = fPhi_default;
  fHitMrdChankeys.clear();
  fCoordsAtLayer.clear();
  fPaddleProbs.clear();
  mrddigitchankeysthisevent.clear();
  h2_angle->Clear();
  h2_space->Clear();
  
  if(!m_data->Stores.count("ANNIEEvent") ){ 
    Log("MrdLikelihoodTracker tool: No ANNIEEvent store!",v_error,fVerbose); 
    return false;
  }

  int get_ok = 0;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("MCEventNum",fEventNumber);
  if(not get_ok){ Log("MrdLikelihoodTracker tool: Error retrieving MCEventNum, true from ANNIEEvent!",v_error,fVerbose); return false; }
  
  std::vector<std::vector<int>> MrdTimeClusters;
  get_ok = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
  if (not get_ok) { Log("MrdLikelihoodTracker Tool: Error retrieving MrdTimeClusters map from CStore, did you run TimeClustering beforehand?",v_error,fVerbose); return false; }

  int nClust = MrdTimeClusters.size();
//  std::cout << "Number of MrdTimeClusters " << nClust << "\n";
  if(nClust != 1) return true; // TODO maybe add a log message
  std::cout << "Size of MrdTimeClusters[0] " << MrdTimeClusters[0].size() << "\n";

  //get_ok = m_data->CStore.Get("MrdDigitTimes",MrdDigitTimes);
//  if (not get_ok) { Log("EventDisplay Tool: Error retrieving MrdDigitTimes map from CStore, did you run TimeClustering beforehand?",v_error,verbose); return false; }
  get_ok = m_data->CStore.Get("MrdDigitChankeys",mrddigitchankeysthisevent);
  if (not get_ok) { Log("EventDisplay Tool: Error retrieving MrdDigitChankeys, did you run TimeClustering beforehand",v_error,fVerbose); return false;}
  
  for(int digit_value : MrdTimeClusters[0]){
//    std::cout << digit_value << "\n";
    unsigned long chankey = mrddigitchankeysthisevent.at(digit_value);
//    Detector *thedetector = fGeom->ChannelToDetector(chankey);
//    unsigned long detkey = thedetector->GetDetectorID();
//    if (thedetector->GetDetectorElement()!="MRD") continue;
//    double mrdtimes=MrdDigitTimes.at(digit_value);
    fHitMrdChankeys.push_back(chankey);    
//    std::cout << chankey << "\n";
  }
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bool fit = true;
  bool draw = true;
  if(fit){
    //  fMinimizer = new TMinuitMinimizer(ROOT::Minuit::kMigrad, fNPars);
  //  std::cout << fMinimizer << "\n";
    //fMinimizer = std::unique_ptr<ROOT::Math::Minimizer> ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );
    fMinimizer = ROOT::Math::Factory::CreateMinimizer( "Minuit", "Migrad" ) ;
  //  std::cout << fMinimizer << "\n";
    
    DefineFunc();
    fMinimizer->SetFunction(fFunc);
  
    fMinimizer->SetMaxFunctionCalls(10000);
    fMinimizer->SetMaxIterations(10000);
    fMinimizer->SetTolerance( 1e-12 );
    fMinimizer->SetPrintLevel(1);
  
    fMinimizer->SetVariable( 0, "start_x", 0., 0.005);
    fMinimizer->SetVariableLimits( 0, -2., 2.); // TODO refine these limits
    
    fMinimizer->SetVariable( 1, "start_y", 0., 0.005);
    fMinimizer->SetVariableLimits( 1, -2., 2.); // TODO refine these limits
    
    fMinimizer->SetVariable( 2, "theta", 45., 0.1);
    fMinimizer->SetVariableLimits( 2, 0.1, 89.9);
  
    fMinimizer->SetVariable( 3, "phi", 10., 0.5);
    fMinimizer->SetVariableLimits( 3, -360., 720.); 

    DoFit();
    
    delete fMinimizer;
    fMinimizer = NULL;

    if(!fFitStatus){
      // I noticed that if there was no succesful fit we were carrying over fit params from
      // previous event and putting it in the store. Here I am instead putting some dummy params. 
      // TODO this is a bad solution, fix it !!
      fFitVals = {0., 20., 88., 0.};
      m_data->Stores["ANNIEEvent"]->Set("MrdLikelihoodFitVals",fFitVals);
      return true;
    }
      
    if(fFitVals.size() == fNPars){
      std::cout << "~~~ Event num " << fEventNumber << " ~~~ \n";
      std::cout << fFitVals.at(0) << "\n";
      std::cout << fFitVals.at(1) << "\n";
      std::cout << fFitVals.at(2) << "\n";
      std::cout << fFitVals.at(3) << "\n";     
      m_data->Stores["ANNIEEvent"]->Set("MrdLikelihoodFitVals",fFitVals);
    }
    else{
      Log("MrdLikelihoodTracker tool: Number of fit params greater than number of model params. Something is wrong here",v_error,fVerbose);
      return false;
    }

  }//end if(fit)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(draw){

    // draw angle
    fStartX = fFitVals.at(0);
    fStartY = fFitVals.at(1);
    fTheta = fFitVals.at(2);
    fPhi = fFitVals.at(3);
    for(int iBinX=1; iBinX<=nBinsX; iBinX++){
      double binX = h2_angle->GetXaxis()->GetBinCenter(iBinX);
      fTheta = binX;
      for(int iBinY=1; iBinY<=nBinsY; iBinY++){
        double binY = h2_angle->GetYaxis()->GetBinCenter(iBinY);
        fPhi = binY;        

        FillCoordsAtLayer();
        FillPaddleProbs();
        
        double prob = 1.;
        for(int chankey : fAllMrdChankeys){          
          bool isHit = false;
          double val = 1.;
          if(std::find(fHitMrdChankeys.begin(), fHitMrdChankeys.end(), chankey) != fHitMrdChankeys.end()){  
            val = fPaddleProbs.at(chankey); 
            isHit = true;
        }
          else{
            val = ( 1.- fPaddleProbs.at(chankey) );           
          }    
          prob *= val;
        }// end loop over all MRD channels
        double val = -2. * std::log(prob);
        h2_angle->SetBinContent(iBinX, iBinY, val);
      }//end loop over y-bins
    }//end loop over x-bins    

    // draw space
    // TODO can't I just use Likelihood() method here?
    fStartX = fFitVals.at(0);
    fStartY = fFitVals.at(1);
    fTheta = fFitVals.at(2);
    fPhi = fFitVals.at(3);
    for(int iBinX=1; iBinX<=nBinsX; iBinX++){
      double binX = h2_space->GetXaxis()->GetBinCenter(iBinX);
      fStartX = binX;
      for(int iBinY=1; iBinY<=nBinsY; iBinY++){
        double binY = h2_space->GetYaxis()->GetBinCenter(iBinY);
        fStartY = binY;        

        FillCoordsAtLayer();
        FillPaddleProbs();
        
        double prob = 1.;
        for(int chankey : fAllMrdChankeys){          
          bool isHit = false;
          double val = 1.;
          if(std::find(fHitMrdChankeys.begin(), fHitMrdChankeys.end(), chankey) != fHitMrdChankeys.end()){
            val = fPaddleProbs.at(chankey); 
            isHit = true;
          }
          else{
            val = ( 1.- fPaddleProbs.at(chankey) );
          }    
          prob *= val;
        }// end loop over all MRD channels
        double val = -2. * std::log(prob);
        h2_space->SetBinContent(iBinX, iBinY, val);
      }//end loop over y-bins
    }//end loop over x-bins    
    DrawIt();

  }//end if(draw)
    
  return true;
}

//..................................................................................

void MrdLikelihoodTracker::DefineFunc()
{
  fFunc = ROOT::Math::Functor(
    [&](double const *pars) {
      double ret = -2. * std::log( Likelihood(pars) ); // negative so we can _minimize_
//      double ret = -1. *  Likelihood(pars); // negative so we can _minimize_ 
      return ret;
    },
    fNPars);
}

//..................................................................................

double MrdLikelihoodTracker::Likelihood(const double *parVals)
{  
  fStartX = parVals[0];
  fStartY = parVals[1];
  fTheta = parVals[2];
  fPhi = parVals[3];
  FillCoordsAtLayer();
  FillPaddleProbs();
  double prob = 1.;
  for(int chankey : fAllMrdChankeys){

    bool isHit = false;
    double val = 1.;
    if(std::find(fHitMrdChankeys.begin(), fHitMrdChankeys.end(), chankey) != fHitMrdChankeys.end()){
      val = fPaddleProbs.at(chankey); 
      isHit = true;
    }
    else{
      val = ( 1.- fPaddleProbs.at(chankey) );
    }    
    prob *= val;    
  }// end loop over all MRD channels  
  
  return prob;
}

//..................................................................................

void MrdLikelihoodTracker::DoFit()
{    
  fFitStatus = fMinimizer->Minimize();
  if(fFitStatus){
    for(int i=0; i<fNPars; i++){
      double val = fMinimizer->X()[i];
      if(i == 3){
        if(val < 0) val = 360. + val;
        else if (val > 360) val -= 360.;
      }
      fFitVals.push_back(val);
    }// end for(fNPars)
  }// end if(fFitStatus)
}

//..................................................................................

void MrdLikelihoodTracker::FillCoordsAtLayer()
{
  fCoordsAtLayer.clear(); // TODO might be unecessary?
  for(int iLayer=2; iLayer < (fNLayer+2); iLayer++){ // TODO not a fan of this hardcoded manual offset
    fCoordsAtLayer[iLayer] = CalcCoordsAtZ( fZ_midpoints.at(iLayer), fStartX, fStartY, fZ_start, fTheta, fPhi);
  }
  
}

//..................................................................................

void MrdLikelihoodTracker::FillPaddleProbs()
{  
  fPaddleProbs.clear(); // TODO might be unecessary?
  
  for(int chankey : fAllMrdChankeys){

    Paddle *mrdpaddle = (Paddle*) fGeom->GetDetectorPaddle(chankey);
    
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

    int iLayer = mrdpaddle->GetLayer();
    
    std::pair<double, double> coords = fCoordsAtLayer.at(iLayer);
    
    //TODO bad name
    double chi2_x = (coords.first - x_mid)*(coords.first - x_mid) / (x_extent * x_extent);
    double chi2_y = (coords.second - y_mid)*(coords.second - y_mid) / (y_extent * y_extent);


    double prob = ( 1./(2*M_PI*x_extent*y_extent) ) * exp( - (chi2_x/2.) - (chi2_y/2.) );
    if(prob > 1.){
      std::cout << "~~~~ \n";
      std::cout << prob << "\n";
      std::cout << chi2_x << " : (" << coords.first << " - " << x_mid << ") / " << x_extent << "\n";
      std::cout << chi2_y << " : (" << coords.second << " - " << y_mid << ") / " << y_extent << "\n"; 
      abort();
    }
    
    fPaddleProbs[chankey] = prob;
  }// end loop over MRD paddles
  
}


//..................................................................................

void MrdLikelihoodTracker::DrawIt()
{
  TCanvas* c = new TCanvas();
//  c->SetLogz(true);
  c->SetRightMargin(0.13);

  double yLo = 0;
  double yHi = 350;
  
  h2_angle->SetStats(0);
  h2_angle->GetZaxis()->SetRangeUser(yLo, yHi);
  h2_angle->GetZaxis()->SetTitle("-2ln(L)");
  h2_angle->GetXaxis()->SetTitle("#theta (degrees)");
  h2_angle->GetYaxis()->SetTitle("#phi (degrees)");
  h2_angle->GetZaxis()->CenterTitle(true);
  h2_angle->GetYaxis()->CenterTitle(true);
  h2_angle->GetXaxis()->CenterTitle(true);

  int iMinBin = h2_angle->FindBin(fFitVals.at(2), fFitVals.at(3));
  double minVal = h2_angle->GetBinContent(iMinBin);
  
  h2_angle->SetTitle(("Angular likelihood surface at best fit spatial points, " + std::to_string(minVal) ).c_str() );
  TMarker* m_angle = new TMarker(fFitVals.at(2), fFitVals.at(3), 29);  
  m_angle->SetMarkerSize(3);
  m_angle->SetMarkerColor(kRed);
  h2_angle->Draw("colz");
  m_angle->Draw("SAME");
  c->SaveAs( ("./MrdPaddlePlots/surface_"+std::to_string(fEventNumber)+"_angle.png").c_str() );

  h2_space->SetStats(0);
  h2_space->GetZaxis()->SetRangeUser(yLo, yHi);
  h2_space->GetZaxis()->SetTitle("-2ln(L)");
  h2_space->GetXaxis()->SetTitle("Start X (m)");
  h2_space->GetYaxis()->SetTitle("Start Y (m)");
  h2_space->GetZaxis()->CenterTitle(true);
  h2_space->GetYaxis()->CenterTitle(true);
  h2_space->GetXaxis()->CenterTitle(true);
  h2_space->SetTitle("Spatial likelihood surface at best fit angle points");
  TMarker* m_space = new TMarker(fFitVals.at(0), fFitVals.at(1), 29);  
  m_space->SetMarkerSize(3);
  m_space->SetMarkerColor(kRed);
  h2_space->Draw("colz");
  m_space->Draw("SAME");
  c->SaveAs( ("./MrdPaddlePlots/surface_"+std::to_string(fEventNumber)+"_space.png").c_str() );
  
  delete c;
  c = nullptr;

}

//..................................................................................

bool MrdLikelihoodTracker::Finalise()
{

  delete h2_angle;
  h2_angle = nullptr;

  delete h2_space;
  h2_space = nullptr;
  
  return true;
}

//..................................................................................
