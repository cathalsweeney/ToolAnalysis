#ifndef MrdLikelihoodTracker_H
#define MrdLikelihoodTracker_H

#include <string>
#include <iostream>

#include "Tool.h"
//#include "Geometry.h"

#include "TCanvas.h"
#include "TH1.h"

/**
 * \class MrdLikelihoodTracker
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class MrdLikelihoodTracker: public Tool {


 public:

  MrdLikelihoodTracker(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.
  bool FillCellProbs();


 private:

    /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int fVerbose;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;

  
  int fEventNumber;
  double fStartX;
  double fStartY;
  double fStartZ;
  double fTheta;
  double fPhi;
  
  Geometry *geom = nullptr;  

  TH1D* hist = nullptr;
};


#endif
