#include "HCVariables.hh"
#include "Rtypes.h"

#include <iostream>
#include "string.h"
#include "string"
#include <fstream>
#include <sstream>

#include "TString.h"

using namespace std;


HCVariables::HCVariables()
{
  cout << "<HCVariables::HCVariables>: LOADING" << endl;
  TString inputfile = "options.ini";
  LoadFromFile(inputfile);

  // I know that there should be a way to call these functions independently of this one
  CreateFolder();
}

HCVariables::~HCVariables()
{
  cout << "<HCVariables::~HCVariables>: END" << endl;
}


int
HCVariables::LoadFromFile(TString FileName)
{
  // The use of this algorithm is much better in C++ style
  // where each variable can be called defining the class
  // and accessing the methods as Variable-><option>
  // For the moment, it uses global variables.
  // Nevertheless, I will read them once and them move them throught
  // the code as arguments.

  
  ifstream infile;
  string line;
  infile.open(FileName);
  if (infile)
    {
      while (!infile.eof())
	{
	  if (getline(infile,line))
	    {
	   
	      string buf;
	      string wert;
	      TString TSwert;
	      int integerwert;
	      double doublewert;
	      string boolewert;
	      
	      stringstream ss(line);
	      //	      cout << ">" << line << endl;
	      ss >> buf;
	      if (buf == "project")	 
	     	{
		  ss >> TSwert;
		  nameProject = TSwert;
		}
	      
	      if (buf == "Ebeam")	 
	     	{			
	      	  ss >> doublewert;
		  gBeam = doublewert;
		}

	      if (buf == "angle")	 
	     	{			
	      	  ss >> doublewert;
		  gTheta = doublewert;
		}
	      if (buf == "p_min")	 
	     	{			
	      	  ss >> doublewert;
		  Ep_min = doublewert;
		}
	      if (buf == "W_threshold")	 
	     	{			
	      	  ss >> doublewert;
		  WThreshold = doublewert;
		}
	      if (buf == "plotPS")	 
	     	{			
	      	  ss >> boolewert;
		  plotPhaseSpace = ReadBoolean(boolewert);
		}
	      if (buf == "plotSF")	 
	     	{			
	      	  ss >> boolewert;
		  plotStructureF = ReadBoolean(boolewert);
		}
	      if (buf == "plotSFrep")	 
	     	{			
	      	  ss >> wert;
		  plotStructureFrep = wert;
		}
	      if (buf == "plotXSHe3")	 
	     	{			
	      	  ss >> boolewert;
		  plotCrossSecHe3 =  ReadBoolean(boolewert);
		}
	      if (buf == "target")	 
	     	{
		  ss >> TSwert;
		  nucleon = TSwert;
		}

	      if (buf == "IPOL")	 
		{
		  ss >> integerwert;
		  varIPOL =integerwert;
	      	}

	      if (buf == "IA1")
	       	{
		  ss >> integerwert;
		  varIA1 =integerwert;
		}
		
	      if (buf == "AsymChoice")
		{
		  ss >> integerwert;
		  varAsym =integerwert;
		}

	      if (buf == "SFChoice")
		{
		  ss >> integerwert;
		  varSFC = integerwert;	
		}

	      if (buf == "YoniIndex")
		{
		  ss >> integerwert;
		  varYoniIndex = integerwert;	
		}
	      if (buf == "reso")
		{
		  ss >> wert;
		  varreso = wert;	
		}
	    }
	}
    }

  Target_rad();
  // with this function we define the masses and type of target once and for all
  // I prefer to separate the function and call it from here. 

  radTheta = gTheta/rad2deg; 

  return 1;
}


// This function is auxiliary to the Variables function
bool HCVariables:: ReadBoolean(string Value)
{
  cout<<"RB: "<<Value<<endl;
  
  if ((Value=="true")
      || (Value=="True")
      || (Value=="TRUE")
      || (Value=="1"))
    {return true;}
  else
    {return false;}
}


void HCVariables::Target_rad()
{
  // Ideally, we should use SWITCH/CASE
  // but switch doesn't handle strings very well

  // Honestly, this function is a bit useless except to select
  // the case, but the mass is redifined when we want to use He3 in the
  // SelectProcess. Definetely, this could be simplified with a C++ style

  if( nucleon == "neutron")
    {
      massTarget = mn;
      typeTarget = 0;
    }
  else  if( nucleon == "proton")
    {
      massTarget = mp;
      typeTarget = 1;
    }
  else  if( nucleon == "He3")
    {
      massTarget = mHe3;
      typeTarget = 3; // in the improbable case we want D2 in the future
    }
  else
    {
      cout << "Invalid Target" << endl;
      exit(0);
    }
        
  //  cout<<"massTarget: "<<mTar<<endl; 

  return;
  
}

void HCVariables::CreateFolder()
{
  // create the folder where to store the data  
  TString prefix = "./"; 
  WorkFolder = prefix + nameProject;
  cout<< nameProject <<endl;
  gSystem->MakeDirectory(WorkFolder);
  
}




HCVariables *Variables=NULL;
