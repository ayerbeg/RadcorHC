#include "HCSelectProcess.hh"
#include "Rtypes.h"

#include <iostream>
#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TString.h"


#include "HCVariables.hh"

using namespace std;



HCSelectProcess::HCSelectProcess()
{
  cout << "<HCSelectProcess::HCSelectProcess>: LOADING" << endl;
}

HCSelectProcess::~HCSelectProcess()
{
  cout << "<HCSelectProcess::~HCSelectProcess>: END" << endl;
}


void HCSelectProcess::SelectProcess()
{
  cout<<"massTarget: "<<Variables->massTarget<<endl;
  cout<<"typeTarget: "<<Variables->typeTarget<<endl;

  TString FinalFolder;
  
  double mn, mp;

  vector<double> vEs_00;
  vector<double> vEp_00;
  
  // Do I really need the FULL Phase Space? 
  
  switch(Variables->typeTarget)
    {
      // neutron and proton phase-space are similar except for the mass
      // WFolder is the folder created for the project (in options.ini)
      // FinalFolder is the specific folder for a given angle (in options.ini)

    case 0: //neutron
    case 1: //proton
      //      PSpace ->create_space_phase(vEs_00, vEp_00); // To be quite precise, there should be space-phase for
      PSpace ->create_space_phase(); // To be quite precise, there should be space-phase for
      // neutron and for protons separately. I think having just proton is good enough.
      // This can be modified easily.
      
      cout<<"HCSP: "<<PSpace->vEs_PS[12]<<" "<<PSpace->vEs_PS[60] <<endl;

      Spectra -> CreateSpectra();
      //    clas_model();

      break;

    case 3:
      // first the neutron
      //   create_space_phase(WFolder, mTarget);
      //   CreateSpectra(WFolder, "neutron", Variables->Ep_min, vEsneu, vEpneu, FinalFolder, vvEsneu, vvEpneu, vvnuneu);

      //    clas_model(FinalFolder, "neutron", mTarget, vvEsneu, vvnuneu, vvEpneu, vvF1neu, vvF2neu);
      
      // second the proton
      //    create_space_phase(WFolder, mTarget);
      //   CreateSpectra(WFolder, "proton", Variables->Ep_min, vEspro, vEppro, FinalFolder, vvEspro, vvEppro, vvnupro);

      //     clas_model(FinalFolder, "proton", mTarget, vvEspro, vvnupro, vvEppro, vvF1pro, vvF2pro);

      // combine to create the Cross-Section
      // It is the same conundrum... I am making the cross-section for all the Es bins
      // in the created phase-space, BUT we are asuming that the bins in proton and neutron
      // are the same, here they are the same, but is it true for every kinematics?
      // I will use the proton bins
     
      //      MakeCross(vvEspro, vvnupro, vvF1neu, vvF2neu, vvF1pro, vvF2pro, WFolder);
      break;
    default:
      
      break;
    }
  
  
}
