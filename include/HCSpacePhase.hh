#ifndef HCSpacePhase_h
#define  HCSpacePhase_h 1

#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TString.h"

#include "HCVariables.hh"

using namespace std;

class HCSpacePhase
{
public:
  HCSpacePhase();
  ~HCSpacePhase();

  
  void create_space_phase();
  void create_space_phase_proton();
  double sinsq();
  double Es_min();
  double Ep_max();
  double Ep(double);

  
 
  double W2threshold;
 
  double Es_max;
  double mNuc;
  int Esbins;

 
  vector<double> vEs_PS;
  vector<double> vEp_PS;

  TString PSinfo;

    
};

extern HCSpacePhase *PSpace; // I should rename ALL the stuff as PhaseSpace

#endif

