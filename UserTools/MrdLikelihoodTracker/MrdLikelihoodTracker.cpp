#include "MrdLikelihoodTracker.h"

MrdLikelihoodTracker::MrdLikelihoodTracker():Tool(){}


bool MrdLikelihoodTracker::Initialise(std::string configfile, DataModel &data){

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
  
  return true;
}


bool MrdLikelihoodTracker::Execute(){
  Log("MrdLikelihoodTracker tool: Executing.",v_message,fVerbose);

  if(!m_data->Stores.count("ANNIEEvent") ){ 
    Log("MrdLikelihoodTracker tool: No ANNIEEvent store!",v_error,fVerbose); 
    return false;
  }

  int this_event;
  int get_ok = 0;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("MCEventNum",this_event);
  if(not get_ok){ Log("MrdLikelihoodTracker tool: Error retrieving MCEventNum, true from ANNIEEvent!",v_error,fVerbose); return false; }
  std::cout << "this_event is " << this_event << "\n";

  std::vector<std::vector<int>> MrdTimeClusters;
  get_ok = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
  if (not get_ok) { Log("MrdLikelihoodTracker Tool: Error retrieving MrdTimeClusters map from CStore, did you run TimeClustering beforehand?",v_error,fVerbose); return false; }

  int nClust = MrdTimeClusters.size();
  std::cout << "Number of MrdTimeClusters " << nClust << "\n";
  if(nClust > 1) return false; // TODO maybe add a log message
  else if(nClust == 1) std::cout << "Size of MrdTimeClusters[0] " << MrdTimeClusters[0].size() << "\n";
  
  return true;
}


bool MrdLikelihoodTracker::Finalise(){

  return true;
}
