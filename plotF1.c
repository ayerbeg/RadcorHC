{
  Double_t rad2deg = 180./(4.*atan(1.));
  Double_t Thetai = 9.6/rad2deg; // min acceptance angle Murchana thesis
  Double_t Theta = Thetai; // angle
  
  string line;
  string it1, it2, it3, it4, it5, it6, it7, it8; // MY FILES HAVE Es, Ep and nu = (Es-Ep)
  Double_t Es, nu, F1, W2, Q2;
  Int_t file_id = 0;
  bool fileexist = true; 
  vector<Double_t> vEs[100];
  vector<Double_t> vnu[100];
  vector<Double_t> vW2[100];
  vector<Double_t> vQ2[100];
  vector<Double_t> vF1[100];
    
  
  // THIS LOOP KEEPS READING THE FILES IN A GIVEN FOLDER (THE ONE DEFDISD BY Ep)
  while(fileexist)
    {
      ifstream  ReadFileDIS(Form("/home/ayerbe/Analysis/RC-Phys-Models/d2n-Carlos/KINA_DIS/%2.2fdeg/XS_%d.dat",Theta *rad2deg, file_id ));
      if (ReadFileDIS)
	{
	  while (!ReadFileDIS.eof())
	    {
	      while ( getline( ReadFileDIS, line ) )
		{
		  stringstream ss( line );      // Set up up a stream from this line
	     
		  ss>>it1>>it2>>it3>>it4>>it5>>it6>>it7>>it8;
		  Es =stod( it1);
		  nu =stod( it2);
		  W2 =stod( it3);
		  Q2 =stod( it4);
		  F1 =stod( it5);

		  vEs[file_id].push_back(Es);
		  vnu[file_id].push_back(nu);
		  vW2[file_id].push_back(W2);
		  vQ2[file_id].push_back(Q2);
		  vF1[file_id].push_back(F1);

		  //	  cout<<"file: "<<file_id<<" Es: "<<Es<<endl;
		}
	    }
	  file_id++;
	}
      else
	{
	  fileexist = false;
	}
    }


  TCanvas *c3;
  c3 = new TCanvas("c3","A Simple Graph Example",200,10,700,500);
  c3->Divide (4,4);
  
  TGraph *grF1[100];
  TMultiGraph *mg = new TMultiGraph();
  
  for (Int_t i=0;i<100;i++) 
    {
      grF1[i] = new TGraph(vnu[i].size(), &vnu[i][0], &vF1[i][0]);	 

    }
  TAxis *axisINE[16]; 
  for (Int_t i=1;i<=16;i++) 
    {
      c3->cd(i);
      grF1[i]  ->SetLineWidth(4);
      grF1[i] ->SetMarkerColor(4);
      grF1[i] ->SetLineColor(kPink);

      axisINE[i] =grF1[i] ->GetXaxis();
      axisINE[i] -> SetLimits(1.5,4.); 
      grF1[i] -> GetHistogram()->SetMaximum(1.2);   // along          
      grF1[i] -> GetHistogram()->SetMinimum(.4);
      
      grF1[i]-> GetHistogram()->SetXTitle("#nu (GeV)");
      grF1[i]-> GetHistogram()->SetYTitle("neutron F1");
      grF1[i]-> Draw("AL");
    }







  
}
