#ifndef HCVariables_h
#define HCVariables_h 1

#include "string.h"
#include "string"
#include <fstream>
#include <sstream>
#include "TSystem.h"
#include "TString.h"

using namespace std;

class HCVariables
{
public:
  HCVariables();
  ~HCVariables();

  int LoadFromFile(TString);
  bool ReadBoolean(string);
  //  TString inputfile;

  void Target_rad();
  void CreateFolder();
  
  const Double_t rad2deg = 180./(4.*atan(1.));
  TString nameProject;
  TString nucleon;
  double gTheta;
  double radTheta;
  double Ep_min;
  double gBeam;
  double WThreshold;
  bool plotPhaseSpace;
  bool plotStructureF;
  bool plotCrossSecHe3;
  string plotStructureFrep;

  
  // var- prefix to work in the code internally
  int varSFC ;
  int varAsym;
  int varIA1;
  int varIPOL; 
  int varYoniIndex;
  string varreso;


  // to the constants (in GeV)
  double  mp = 0.938272 ;//proton mass
  double  mn = 0.939565 ;//neutron mass
  double  mHe3 = 2.81 ;// GeV. 3 uma
  double  mpion = 0.140 ;// pion mass 

  double massTarget;
  int typeTarget;



  TString WorkFolder; //So I can access to the workfolder anytime
  
};

extern HCVariables *Variables;

#endif
