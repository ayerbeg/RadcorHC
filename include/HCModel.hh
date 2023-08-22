#ifndef HCModel_h
#define  HCModel_h 1

#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TString.h"

#include "HCVariables.hh"
#include "HCSpacePhase.hh"
#include "HCSpectra.hh"
using namespace std;


extern"C"{
  void strucfunc_(int *, double *, double *, char *, int *, int *, int *, int *, int *, char *, double *, double *, double *, double *, double *, double *);
}

class HCModel
{
public:
  HCModel();
  ~HCModel();

// Order and declaration of the arguments:
  // (i)file_idx, (d)W2, (d)Q2, (c)Target, (i)YoniIndex,
  // (i)IPOL, (i)IA1, (i)SFC, (i)Asym, (c)reso,
  // (d)F1, (d)F2, (d)g1, (d)g2, (d)A1, (d)A2
  

  void model2use();
  void clas_model(TString);
  void MakeCrossHe3();// to generalize for other targets, maybe could be useful use it with a parameter
  double He3CrossSection( Double_t,  Double_t, Double_t, Double_t,  Double_t);
  double calcQ2(double, double, double);
  double calcW2(double, double, double);
  



  double Theta;
 

  Double_t Es, Ep, nu;


  double targetMass;
  // I am not sure about these variables
  // I declare them to keep the methods separete
  // With a bit better coding, a lot of  loops could be avoinded
  Double_t Estmp, Eptmp, nutmp;

  // SF data/spectra
  vector<double> vF1;
  vector<double> vF2;

  // THE RETURNED VALUES ACCORDING TO DARREN UPTON DOCUMENT
  // MAYBE WE DON'T NEED ALL OF THEM FOR THE RADIATIVE CORRECTIONS
  // BUT WE CAN KEEP THEM FOR FURTHER ANALYSIS

  Double_t Q2;  // transfer momentum
  Double_t W2;  // invariant mass
  Double_t xbj; // x Bjorken
  Double_t F1=0;  // unpolarized SF1
  Double_t F2=0;  // unpolarized SF2
  Double_t R;   // ?
  Double_t A1;  // asymmetry A1
  Double_t A2;  // asymmetry A2 (?)
  Double_t g1=0;  // polarized SF1
  Double_t g2=0;  // polarized SF2


  
  vector<vector<double>> vvF1pro;
  vector<vector<double>> vvF2pro;

  vector<vector<double>> vvF1neu;
  vector<vector<double>> vvF2neu;


  //The Make Cross section
  
  //  Double_t Es, nu, Ep,  W2, Q2, g1, g2; // Es and nu are common for n and p, W2 and Q2 are not used. g1 and g2 are in stand-by
  Double_t pF1, pF2, nF1, nF2;
  
  Double_t XS;
  
  vector<double> vXS;
  vector<vector<double>> vvXS;

  
};

extern HCModel *Model; // I should rename ALL the stuff as PhaseSpace

#endif


  
