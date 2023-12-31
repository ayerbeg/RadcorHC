#ifndef HCSelectProcess_h
#define HCSelectProcess_h 1

#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TString.h"
#include "HCSpacePhase.hh"
#include "HCSpectra.hh"

using namespace std;

class HCSelectProcess
{
public:
  HCSelectProcess();
  ~HCSelectProcess();

  void SelectProcess();

    
  vector<double> vEs;
  vector<double> vEp;

  vector<double> vEspro;
  vector<double> vEppro;

  vector<double> vEsneu;
  vector<double> vEpneu;


  vector<vector<double>> vvEs;
  vector<vector<double>> vvEp;
  vector<vector<double>> vvnu;

  // I believe that these vectors (pro/neu)
  // content the same or similar values. 
  // (update) Actuallym, they have similar values due to the
  // round to the 2nd decimal. For sanity, I'll keep them separate.
  vector<vector<double>> vvEspro;
  vector<vector<double>> vvEppro;
  vector<vector<double>> vvnupro;
  

  vector<vector<double>> vvEsneu;
  vector<vector<double>> vvEpneu;
  vector<vector<double>> vvnuneu;

  // The SF spectra 
  vector<vector<double>> vvF1;
  vector<vector<double>> vvF2;

  vector<vector<double>> vvF1neu;
  vector<vector<double>> vvF2neu;

  vector<vector<double>> vvF1pro;
  vector<vector<double>> vvF2pro;


};

#endif
