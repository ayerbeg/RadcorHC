#include "HCSpectra.hh"
#include "Rtypes.h"

#include <iostream>
#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TString.h"
#include "HCSpacePhase.hh"
#include "HCSelectProcess.hh"

using namespace std;



HCSpectra::HCSpectra()
{
  cout << "<HCSpectra::HCSpectra>: LOADING" << endl;
  // PSpace  = new HCSpacePhase() ; 

}

HCSpectra::~HCSpectra()
{
  cout << "<HCSpectra::~HCSpectra>: END" << endl;
}


void HCSpectra::CreateSpectra()
{
  cout<<"*******create spectra*************"<<endl;
 
  cout<<"Theta: "<<Variables->gTheta<<endl;


  // ALL THESE FOLDERS SHOULD BE MOVE TO A DIFFERENT METHOD
   
  //  TString AngFolder = Variables->WorkFolder+ angle_fol+"deg/";
  AngFolder = Variables->WorkFolder+ Form("/%2.2fdeg/", Variables->gTheta);
  NucFolder = AngFolder+Variables->nucleon;

  cout<<"AngFolder"<<AngFolder<<endl;
  cout<<"NucFolder"<<NucFolder<<endl;
  
  gSystem->MakeDirectory(AngFolder);
  gSystem->MakeDirectory(NucFolder); 

  // The NucFolder is the FinalFolder used internally in the rest of the code
  // it is specific for the angle set in options.ini
  // Perhaps it should be improved for similar angle but other characteristics
  // as is right now, it will overwrite any data present.

  // Keeping the philosophy of working without reading any of the created files
  // we need to store the output here in containers to send to the next method
  // but we save the data as was done with the single scripts

   
  Int_t file_id = 0;
  Int_t file_idx =0;

  
  // These vectors are the limits of the space-phase
  // the Es binning and the corresponding Ep according to
  // the W threshold.

  NoEsBins = PSpace->vEs_PS.size();
  Double_t delta = 0;

  cout<<NoEsBins<<" "<<endl;
  
  gSystem->MakeDirectory(NucFolder+"/spectra");
  
  for (int idx = 0; idx< NoEsBins; idx++)
    {
      Es = PSpace->vEs_PS[idx];
      //  cout<<vEs[idx];
      Ep = PSpace->vEp_PS[idx];
      EpDelta = Ep;
      
      fileDIS = NucFolder+"/spectra"+Form("/%d.dat", file_id );
      spectra.open(fileDIS);
      
      while( EpDelta >=(Variables->Ep_min -(0.001) ) )
	{
	  nu = Es - EpDelta;
	  vEs.push_back(Es);
	  vEp.push_back(EpDelta);
	  vnu.push_back(nu);
	  
	  spectra<< Es<<" "<<EpDelta<<" "<<nu<<"\n";
	  delta++;
	  EpDelta = Ep - 0.01* delta;
	}

      // There are two ways to pass the data to the clas model
      // I can pass it as the spectrum is created
      // so, just need to call the function here and pass
      // the vector with the info of the present spectra
      // The other option, used here, is to pass the whole information
      // Probably, in memory usage, the first is more efficient
      
      vvnuCS.push_back(vnu);
      vvEsCS.push_back(vEs);
      vvEpCS.push_back(vEp);

      vEs.clear();
      vEp.clear();
      vnu.clear();
      
      spectra.close();
      file_id++;
      // clas_model method should be called here.
      // the arguments should be vectors of Es, Ep(steps), nu, file_id
      delta = 0;

    }

  cout<<"size: "<<vvnuCS.size()<<endl;
  
}

double HCSpectra::Theta()
{
  return Variables->gTheta/rad2deg; 
}

HCSpectra *Spectra=NULL;
