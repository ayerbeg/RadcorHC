#include "HCModel.hh"
#include "Rtypes.h"

#include <iostream>
#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TString.h"

using namespace std;



HCModel::HCModel()
{
  cout << "<HCModel::HCModel>: LOADING" << endl;
 
}

HCModel::~HCModel()
{
  cout << "<HCModel::~HCModel>: END" << endl;
}

void HCModel::model2use()
{
  cout<<"***********model to use***************"<<endl;

 switch(Variables->typeTarget)
    {
      // neutron and proton phase-space are similar except for the mass
     
    case 0: //neutron
    case 1: //proton
   
      clas_model(Variables->nucleon);

      break;

    case 3:
      // WITH SPECTRA, WE REFER TO THE BINNING OF THE PHASE-SPACE.
      // SINCE THE PHASE-SPACE IS THE SAME FOR NEUTRON AND PROTON (is that correct?)
      // THE BINNING IS THE SAME FOR BOTH NUCLEONS.

      clas_model("proton");
      clas_model("neutron");

      MakeCrossHe3();

      // combine to create the Cross-Section
      // It is the same conundrum... I am making the cross-section for all the Es bins
      // in the created phase-space, BUT we are asuming that the bins in proton and neutron
      // are the same, here they are the same, but is it true for every kinematics?
          
      //      MakeCross(vvEspro, vvnupro, vvF1neu, vvF2neu, vvF1pro, vvF2pro, WFolder);
      break;
    default:
      
      break;
    }
  
  
}


void HCModel::clas_model(TString nucleon)
{

  cout<<"***********CLAS model***************"<<endl;
  cout<< Spectra -> NucFolder<<endl;


  //  Theta = Variables->gTheta/rad2deg; // I think I will put this value in Variables, is in every class
  cout<<"Theta: "<<Variables->gTheta<<endl;


  // THESE LINES ARE STRONGLY RELATED WITH THE CLAS SUBROUTINE
  // I SHOULD TAKE CARE MODIFYING THEM

  // These lines are necesary to simplify the conversion from
  // string to char needed by the fortran subroutine
  char creso;
  if( (Variables->varreso == "y")||(Variables->varreso == "Y"))
    {
      creso = 'y' ;
    }
  else
    {
      creso = 'n' ;
    }
  
  char Target; // A lot to understand here.

  int typeTarget;
  //  Target_rad(mNuc, typeTarget); // it doesn't work when the target is He3

  const char* nucle;  
  // Same as before (remember we are using global variables)
  if(nucleon == "proton")
    {
      Target = 'P';
      nucle = "p";
      targetMass = Variables->mp;
    }
  else
    {
      Target = 'N';
      nucle = "n";
      targetMass = Variables->mn;
    }

  
  cout<<"Options: "<<  Variables->varIPOL<<" "<<Variables->varIA1<<" "<< Variables->varAsym<<" "<<Variables->varreso<<" "<<Target<<endl;


  Int_t file_id = 0;

  bool fileexist = true;

  TString SaveData = Form(Spectra -> NucFolder+"/SF%s",nucle);
  gSystem->MakeDirectory(SaveData);

  cout<< SaveData<<endl;  

  int NoOfSpectra =  Spectra-> vvEsCS.size();// somehow, this is the number of spectra (the file id)
   
  for(int spectraIDX = 0; spectraIDX < NoOfSpectra; spectraIDX++)
    {
   
      int NoOfPoints = Spectra-> vvEsCS[spectraIDX].size();// this is the number of data points in the spectra
      file_id = spectraIDX;

      ofstream SaveFile(SaveData+Form("/%sSF_%d.dat",Variables->radTheta,  nucle, file_id ));
 
      for(int SpecPointsIDX = 0; SpecPointsIDX < NoOfPoints; SpecPointsIDX++)
	{

	  Es = Spectra-> vvEsCS[spectraIDX][SpecPointsIDX];
	  Ep = Spectra-> vvEpCS[spectraIDX][SpecPointsIDX]; 
	  nu = Spectra-> vvnuCS[spectraIDX][SpecPointsIDX]; 
  

	  Q2 = calcQ2(Es,Ep, Variables->radTheta );
	  W2 = calcW2(targetMass, nu, Q2);

	  
	  strucfunc_(&file_id, &W2, &Q2, &Target, &(Variables->varYoniIndex), &(Variables->varIPOL), &(Variables->varIA1), &(Variables->varSFC), &(Variables->varAsym), &creso, &F1, &F2, &g1, &g2, &A1, &A2 );

	  // with the option 'fixed' every value carries the same number of digits, filling with zeroes when necessary
	  // Data saved: Es -- nu -- W2 -- Q2 -- F1 -- F2 -- g1 -- g2
	  SaveFile<<fixed<<Es<<"\t"<<nu<<"\t" << W2 <<"\t"<<Q2<<"\t"<< F1<<"\t"<< F2<<"\t"<< g1<<"\t"<< g2 << endl;
	  // cout<<fixed<<Es<<"\t"<<nu<<"\t" << W2 <<"\t"<<Q2<<"\t"<< F1<<"\t"<< F2<<"\t"<< g1<<"\t"<< g2 << endl;

	  vF1.push_back(F1);
	  vF2.push_back(F2);
	}

      // I overusing the name of the variables, nucleon is a local variable here
    if(nucleon == "proton")
      {
	vvF1pro.push_back(vF1);
	vvF2pro.push_back(vF2);
      }
    else
      {
	vvF1neu.push_back(vF1);
	vvF2neu.push_back(vF2);
      }
	
      vF1.clear();
      vF2.clear();
      SaveFile.close();
     }

 
  // can we combine the previous work in order to calculate THIS:
  //  XS = CrossSection(Es, nu, Theta, 2*pF1+nF1,  2*pF2+nF2);
  // probably it needs an extra function where it call the method twice

}





//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-



//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

// suffix MC for MakeCross (this function)
void HCModel::MakeCrossHe3()
{

  cout<<"*********** Make Cross ***************"<<endl;

  TString SaveData = (Spectra -> NucFolder+"/He3CrossSection");

  gSystem->MakeDirectory(SaveData);
  cout<<SaveData<<endl;  

  int NoOfSpectra = Spectra -> vvEsCS.size();// somehow, this is the number of spectra (the file id)
  cout<<"Hello "<< NoOfSpectra<<endl;

  cout<<  Spectra ->vvnuCS.size()<<" "<<vvF1neu.size()<<" "<<vvF2neu.size()<<" "<<vvF1pro.size()<<" "<<vvF1pro.size()<<endl;
  
  Int_t file_id = 0;
  for(int spectraIDX = 0; spectraIDX < NoOfSpectra; spectraIDX++)
    {
 
      int NoOfPoints =  Spectra ->vvEsCS[spectraIDX].size();// this is the number of data points in the spectra
      file_id = spectraIDX;

      ofstream  SaveFile(SaveData +Form("/XS%2.2f_%d.dat", Variables->gTheta, file_id ));
 
      for(int SpecPointsIDX = 0; SpecPointsIDX < NoOfPoints; SpecPointsIDX++)
	{

	  // All units are in GeV, but the cross-section is handled in MeV
	  // so, we convert all to MeV (1e3) but at the end of the process

	  Estmp = Spectra-> vvEsCS[spectraIDX][SpecPointsIDX];
	  Eptmp = Spectra-> vvEpCS[spectraIDX][SpecPointsIDX]; 
	  nutmp = Spectra-> vvnuCS[spectraIDX][SpecPointsIDX]; 
	  
	  // Es = vvEsMC[spectraIDX][SpecPointsIDX]; 
	  // nu = vvnuMC[spectraIDX][SpecPointsIDX]; 

	  
	  Q2 = calcQ2(Estmp, Eptmp, Variables->radTheta );
	  //What should be the mass? the nucleon? the He3?
	  // I use the mass of the proton
	  W2 = calcW2(Variables->mp, nutmp, Q2)*pow(1e3,2);  //in MeV^2
	  
	  pF1 = vvF1pro[spectraIDX][SpecPointsIDX];   
	  pF2 = vvF2pro[spectraIDX][SpecPointsIDX]; 
	  
	  nF1 = vvF1neu[spectraIDX][SpecPointsIDX];  
	  nF2 = vvF2neu[spectraIDX][SpecPointsIDX];  
	  
	  XS = He3CrossSection(Es*1e3, nu*1e3, Variables->radTheta, 2*pF1+nF1,  2*pF2+nF2);

	  // Table --> Es ,nu, W, XS, 0
	  SaveFile<< Estmp*1e3 << "\t" <<nutmp*1e3  << "\t"<< sqrt(W2) << "\t" <<XS<<"\t"<<0.0<<endl;
	  //	  cout<<fixed<< Es*1e3 << "\t" <<nu*1e3  << "\t"<< sqrt(W2) << "\t" <<XS<<"\t"<<0.0<<endl;
	  vXS.push_back(XS);
	}
      vvXS.push_back(vXS);
      vXS.clear();
      SaveFile.close();
    }

  
  // cout<<"bool: "<<plotCrossSecHe3<<endl;
  // if( plotCrossSecHe3)
  //   PlotHe3XS(vvnuMC[NoOfSpectra-1], vvXS[NoOfSpectra-1]);

}


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

double HCModel::He3CrossSection( Double_t tE,  Double_t tnu, Double_t  tAngle, Double_t  tF1,  Double_t tF2)
{
  // tE incoming energy, tnu energy difference, tAngle scattering angle, tF1 = F1n + 2F1p (same F2)
  // units converted already in MeV from GeV


  // I think I missed a factor cos^2 in the formula
  
  double alpha = 1/137.;
  double hbarc = 197.32698;
  double mott = pow(alpha, 2)/(4*(tE*tE)*pow(sin(tAngle/2.), 4));
  double cross =pow(cos(tAngle/2.), 2)* ((tF2/(tnu))+((2*tF1*pow(tan(tAngle/2.), 2))/Variables->mp));
  
  return mott*cross*(hbarc*hbarc)*1e7;
  // the convertion factor from 1/MeV^3 -> nb/MeV
}


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

double HCModel::calcQ2(double Ei, double Ef, double Th)
{
  return Ei*Ef*4*sin(Th/2.)*sin(Th/2.);
}
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
double HCModel::calcW2(double mp, double nulocal, double Q2local)
{
  return mp*mp + (2*mp*nulocal) - Q2local;
}
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-








HCModel *Model=NULL;
