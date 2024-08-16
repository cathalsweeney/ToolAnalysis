#ifndef MrdLikelihoodTracker_H
#define MrdLikelihoodTracker_H

#include <string>
#include <iostream>

#include "Tool.h"
//#include "Geometry.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
//#include "Minuit2/Minuit2Minimizer.h"
//#include "TMinuitMinimizer.h"


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
  
  
 private:

  void FillPaddleProbs(); ///< Fill fPaddleProbs
  void FillCoordsAtZ(); ///< Fill fCoordsAtZ
  void DoFit();
  void DefineFunc();
  double Likelihood(const double *parVals);

  static constexpr int fNPars = 4;
//  std::unique_ptr<ROOT::Math::Minimizer> fMinimizer = NULL;
  ROOT::Math::Minimizer* fMinimizer = NULL;
//  ROOT::Minuit2::Minuit2Minimizer*  fMinimizer = NULL;
//    ROOT::Minuit::MinuitMinimizer*  fMinimizer = NULL;
//  ROOT::Minuit::TMinuitMinimizer*  fMinimizer = NULL;

  ROOT::Math::Functor fFunc;

  int fFitStatus = 0;
  double fFitVals[fNPars];
  
    /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int fVerbose;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;

  
  int fEventNumber;
  double fStartX;
  double fStartY;
//  double fStartZ;
  double fTheta;
  double fPhi;

  std::vector<int> fHitMrdChankeys;
  
  Geometry *fGeom = nullptr;  

//  TH1D* hist = nullptr;

  std::vector<double> fZ_midpoints=
  { 3.3638, 3.4884, 3.613, 3.7376,
    3.8622, 3.9833, 4.1044, 4.2255,
    4.3466, 4.4677, 4.5888}; ///< z-coordinates of the midpoints of each plane

  double fZ_start = 3.3608; ///< z-coordinate of the start of the MRD

  std::map<double, std::pair<double,double> > fCoordsAtZ; ///< (x,y) coords of track at z midpoint of each plane given current track params; fCoordsAtZ[z_val] = {x,y} 
  std::map<int, double> fPaddleProbs; ///< Probability of this paddle being hit given current track params; fPaddleProbs[chankey] = prob

  std::vector<unsigned long> mrddigitchankeysthisevent;  //from TimeClustering tool
  TH2D* h2_angle;
  TH2D* h2_space;
  int nBinsX = 100;
  int nBinsY = 100;

  std::vector<int> fAllMrdChankeys;
  
};


#endif
