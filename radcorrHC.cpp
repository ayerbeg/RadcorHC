// NOTES:
// The output folders must be constructed previously to run the script
#include "Rtypes.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"

#include <TROOT.h>
#include <TStyle.h>
#include <vector>
#include <math.h>
#include "TStopwatch.h"
#include <sstream>
#include "TSystem.h"
#include "TString.h"

#include <iomanip> //library for precision()



#include "HCVariables.hh"
#include "HCSelectProcess.hh"
#include "HCSpacePhase.hh"
#include "HCSpectra.hh"
#include "HCModel.hh"
#include "HCPlot.hh"

using namespace std;



//*******************************************
// GLOBAL VARIABLES (MAYBE MOVE TO A HEADER)

bool plotPhaseSpace = 0;
bool plotStructureF = 0;
bool plotCrossSecHe3 = 0;
//stringVariables-> plotStructureFrep;


//*******************************************



/*
Double_t round_number2(Double_t);
Double_t round_number3(Double_t);


void Plot_phase_space(vector<Double_t>, vector<Double_t>, double, double, double, double, double);

void variable(TString );
bool ReadBoolean(string);
void Target_rad(double &, int &);

//void SelectProcess(double, int, TString);


*/

const Double_t rad2deg = 180./(4.*atan(1.));

// NOTEL start to remove the round functions. I am not sure how they will carry the further calculations
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-







//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-



  
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


int main(int argc, char *argv[])
{
 
  TApplication *theApp = new TApplication("app", &argc, argv); //<======

  cout << "Compiled on: " __DATE__ " " __TIME__ "." << endl;

  TStopwatch timer;
  timer.Start();

  Variables = new HCVariables();
  PSpace = new  HCSpacePhase();
  Spectra = new HCSpectra();
  Model = new HCModel();
  
  
  // TString inifile = "options.ini";

  gStyle->SetOptFit(1);
  gStyle->SetStatFormat("6.6g");
  
  PSpace ->create_space_phase();
   
  Spectra -> CreateSpectra();
  Model -> model2use();

  if(Variables -> plotPhaseSpace ||
     Variables -> plotStructureF ||
     Variables -> plotCrossSecHe3||
     Variables -> plotHe3XS3D )
    {
      Plot = new HCPlot();
      Plot-> SelectPlot();
    }

  
  timer.Stop();
  
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();

  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

  cout << "-----> CTRL+C to end" << endl;
  
  theApp->Run();

  delete Variables;
  
  return 1;
}







//*********************************************************************************************
//                     AUXILIARY FUNCTIONS (NOT IMPORTANT FOR THE ANALYSIS)
//*********************************************************************************************

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
// In principle the numbers rounded can be stored directly in the file
// but they wwill be needed in a calculus later on.

Double_t round_number2(Double_t Num)
{
  // setprecision() is a great tool to round the number up/down
  // but it only works for stream (like cout). With some C++ -fu
  // the output is forward to a string (stringstream) and then converted
  // to double precision with 2 significant figures (in this case)
  // I choose 2 figures, because the steps to calculate nu=E-E' will
  // be in 0.1GeV steps

  stringstream ss;

  ss<< setprecision(2) << Num<< endl; 
  Double_t rounded = stod(ss.str());

  cout<<"ss: "<<rounded<<endl;    

  return rounded;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

Double_t round_number3(Double_t Num)
{
  // same as the previous function but to three significant figures
  // since the steps in Es for a given Ep bin will be in 

  stringstream ss;

  ss<< setprecision(3) << Num<< endl; 
  Double_t rounded = stod(ss.str());

  // cout<<"ss: "<<rounded<<endl;    

  return rounded;
}
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
