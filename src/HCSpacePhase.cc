#include "HCSpacePhase.hh"
#include "Rtypes.h"

#include <iostream>
#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TString.h"

using namespace std;



HCSpacePhase::HCSpacePhase()
{
  cout << "<HCSpacePhase::HCSpacePhase>: LOADING" << endl;
  W2threshold = pow(Variables->WThreshold,2); 
  Es_max = Variables->  gBeam ;

  // the phase-space folder
  PSinfo = Variables->WorkFolder+"/PS_limits";
  gSystem->MakeDirectory(PSinfo);   
}

HCSpacePhase::~HCSpacePhase()
{
  cout << "<HCSpacePhase::~HCSpacePhase>: END" << endl;
}

void HCSpacePhase::create_space_phase()
{
  create_space_phase_proton();// nothing, just call the function I renamed
}

// I gave the name proton expecting to do proton-neutron
// phase-speace differently. But I think it is good enough having
// just the proton one
void HCSpacePhase::create_space_phase_proton()
{
  cout<<"********create space phase*****************"<<endl;

  cout<<"Theta: "<<Variables->gTheta<<endl;

  mNuc = Variables-> mp;  

  // the phase-space data
  ofstream outfile;
  outfile.open(Form(PSinfo+"/Theta.%.2f.Eprime.%.2f.dat", Variables->gTheta, Variables->Ep_min  ) );  
 
  cout<<"mNuc: "<<mNuc<<endl;
  
  // My intention is to create a general code, but it is greatly dependent of certain values
  // like the energy beam. For instance we want 0.1 GeV space bins in Es, but the last bin, can only
  // 0.08 so it needed to be hard coded. Maybe, we can make a code to determine when the division is
  // not exact and add the last bin.
  
  // we need to calculate the number of bins from Es_min and Es_max (10.300) with 0.1 GeV spacing
  // because we round Es_min up to two significatn figures, we assure that the number of bins
  // is integer, not just because we defined it integer (that would truncate the number)
  // Last bin, as suggested by Melanie's theses, will be of 0.8 and added later

  double EsBinSize = 0.1;
  Esbins = ((Es_max - 0.080 ) - Es_min()  ) /(EsBinSize );


  // the number of bins here, must add 0 index and the extra bin at the end
  // due to the irregular bin of the Beam energy
  cout<<"HCSP::ESbins: "<<Esbins<<endl;
  
  //*****************************************************************
  // These lines create the spectrum of Ep vs Es in steps of 0.1 GeV
  // If we want two decimals precisiom in Es, we should add a round
  // function
  //*****************************************************************
  for (int Esidx = 0; Esidx <= Esbins; Esidx++)
    {
      vEs_PS.push_back(Es_min() + EsBinSize *Esidx);
      
      vEp_PS.push_back( Ep(vEs_PS[Esidx] ) );
    }
  
  // This extra line is to include the last bin which won't fit the bin separation criteria odf 0.1GeV
  vEs_PS.push_back(Es_max);
  vEp_PS.push_back( Ep(vEs_PS[Esbins+1])) ;
  
  	  
  //*********************************************************************************
  // These lines store in a ASCII file the Ep-bin vs Es-bin for each Variables->Ep_min and Theta
  // Also creates the folder with the .dat files with the steps to feed the models
  for (int gr_idx0 = 0; gr_idx0 < vEs_PS.size(); gr_idx0++)
    {
      outfile << fixed << vEs_PS[gr_idx0] <<" " << vEp_PS[gr_idx0] <<"\n";
    }
  //*********************************************************************************

  // cout<<"Es size: "<< vEs_PS.size()<<endl;  
    outfile.close();

  return;
}


double HCSpacePhase::Es_min()
{
return (W2threshold - pow(mNuc,2) + 2*mNuc*Variables->Ep_min )/( 2*mNuc - 4*Variables->Ep_min  *sinsq() ) ;
}

double HCSpacePhase::Ep_max()
{
  return ( pow(mNuc,2) + 2*mNuc*Es_max - W2threshold )/( 2*mNuc + 4*Es_max *sinsq() );
}

double HCSpacePhase::Ep(double Es_idx)
{
  return ( pow(mNuc,2) + 2*mNuc*Es_idx - W2threshold)/( 2*mNuc + 4*Es_idx*sinsq() );
}

double HCSpacePhase::sinsq()
{
  return pow(sin((Variables->radTheta )/2.),2); //sin^2(theta/2) --> just to simplify
}



HCSpacePhase *PSpace=NULL;

