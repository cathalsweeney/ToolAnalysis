#include "MrdLikelihoodTracker.h"


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
  if(!m_variables.Get("StartZ",fStartZ)){
    Log("MrdLikelihoodTracker tool: Did not find \"StartZ\" in config file. Quitting",v_message,fVerbose); 
    return false;
  }

  
  if(!m_variables.Get("Theta",fTheta)){
    Log("MrdLikelihoodTracker tool: Did not find \"Theta\" in config file. Quitting",v_message,fVerbose); 
    return false;
  }  
  if(!m_variables.Get("Phi",fPhi)){
    Log("MrdLikelihoodTracker tool: Did not find \"Phi\" in config file. Quitting",v_message,fVerbose); 
    return false;
  }


  hist = new TH1D("hist", "hist", 200, 3., 5.);
  
  
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
  
//  std::vector<std::vector<int>> MrdTimeClusters;
//  get_ok = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
//  if (not get_ok) { Log("MrdLikelihoodTracker Tool: Error retrieving MrdTimeClusters map from CStore, did you run TimeClustering beforehand?",v_error,fVerbose); return false; }
//
//  int nClust = MrdTimeClusters.size();
//  std::cout << "Number of MrdTimeClusters " << nClust << "\n";
//  if(nClust > 1) return false; // TODO maybe add a log message
//  else if(nClust == 1) std::cout << "Size of MrdTimeClusters[0] " << MrdTimeClusters[0].size() << "\n";

  FillCellProbs();

  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
  int n_mrd_pmts = geom->GetNumDetectorsInSet("MRD");
  std::cout << "n_mrd_pmts " << n_mrd_pmts << "\n";
  std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = geom->GetDetectors();

  for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("MRD").begin();
                                                    it != Detectors->at("MRD").end();
                                                  ++it){
    Detector* amrdpmt = it->second;
    unsigned long detkey = it->first;
    unsigned long chankey = amrdpmt->GetChannels()->begin()->first;
    Paddle *mrdpaddle = (Paddle*) geom->GetDetectorPaddle(detkey);
    
    double x_min = mrdpaddle->GetXmin();
    double x_max = mrdpaddle->GetXmax();
    double y_min = mrdpaddle->GetYmin();
    double y_max = mrdpaddle->GetYmax();
    double z_min = mrdpaddle->GetZmin();
    double z_max = mrdpaddle->GetZmax();
    int orientation = mrdpaddle->GetOrientation();    //0 is horizontal, 1 is vertical
    int half = mrdpaddle->GetHalf();                  //0 or 1
    int side = mrdpaddle->GetSide();

    double z_mid = (z_min + z_max)/2.0;
    std::cout << z_mid << "\n";
    hist->Fill(z_mid);
  }  
  
  return true;
}

bool MrdLikelihoodTracker::FillCellProbs()
{
//  double x_0 = 2.;
//  double y_0 = 3.5;
//  double z_0 = 0.;
//  
//  double theta = 10.;
//  double phi = 45.;


  std::vector<double> vCellCenters_z = {0., 2., 4., 6.};

  double sin_phi = std::sin(fPhi*(M_PI/180.) );
  double cos_phi = std::cos(fPhi*(M_PI/180.) );
  double tan_theta = std::tan(fTheta*(M_PI/180.) );
  
  for(double z_val : vCellCenters_z){
    double x_val = fStartX + (z_val - fStartZ)*tan_theta*cos_phi;
    double y_val = fStartY + (z_val - fStartZ)*tan_theta*sin_phi;
    std::cout << "{ " << x_val << ", " << y_val << ", " << z_val << "} \n";
  }
  return true;
}


bool MrdLikelihoodTracker::Finalise()
{
  TCanvas* c = new TCanvas();
  hist->Draw("hist");
  c->SaveAs("foo.png");

  delete hist;
  hist = nullptr;

  delete c;
  c = nullptr;
  
  return true;
}
