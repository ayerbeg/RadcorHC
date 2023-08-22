#ifndef HCPlot_h
#define  HCPlot_h 1

#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"

#include "HCVariables.hh"
#include "HCSpacePhase.hh"
#include "HCModel.hh"
#include "HCSpectra.hh"

using namespace std;

class HCPlot
{
public:
  HCPlot();
  ~HCPlot();

  void SelectPlot();
  
  void Plot_phase_space();

  void PlotSF();
  void GraphSF();
  
  void PlotHe3XS();


  TGraph *grDIS;
  TAxis *axisDIS; 
  TLine *DISlineH;
  TLine *DISlineV;
  
  TLine *INElineStepV[100];
  TLine *INElineStepH[100];

  double Es_max;
  double Es_min;
  double Ep_max;
  double Ep_min;
  vector<Double_t> vgrEsDIS;
  vector<Double_t> vgrEpDIS;

  double eBeam;
  const char* nucle;  
  int lastbin;
  vector<Double_t> vEpPlot;
  vector<Double_t> vnuPlot;
  vector<Double_t> vvF1Plot; //TEMPORARY NAME
  vector<Double_t> vvF2Plot; //TEMPORARY NAME

  TCanvas *c3;
  TCanvas *c4;


  TCanvas *cHe3XS;
  TGraph *grXS;
  vector<double> vXS;
};

extern HCPlot *Plot;

#endif

