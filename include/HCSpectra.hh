#ifndef HCSpectra_h
#define HCSpectra_h 1

#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TString.h"

#include "HCVariables.hh"
#include "HCSpacePhase.hh"

using namespace std;

class HCSpectra
{
public:
  HCSpectra();
  ~HCSpectra();


  void CreateSpectra();
  //HCSpacePhase *sp;
  // HCSpacePhase *PSpace;
  // HCSpacePhase sp; 
  double Theta();
  const Double_t rad2deg = 180./(4.*atan(1.));

  TString AngFolder;
  TString NucFolder;
  
  ofstream spectra;
  double Es, Ep;
  double nu, EpDelta;

  TString fileDIS;

  int NoEsBins;
  
  // These vectors contain the data of a given Es bin
  vector<double> vEs;
  vector<double> vEp;
  vector<double> vnu;


  vector<vector<double>>vvnuCS;
  vector<vector<double>>vvEsCS;
  vector<vector<double>>vvEpCS;
  
  
};
extern HCSpectra *Spectra; // I should rename ALL the stuff as PhaseSpace

#endif
