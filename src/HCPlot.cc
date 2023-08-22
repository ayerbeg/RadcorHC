#include "HCPlot.hh"
#include "Rtypes.h"

#include <iostream>
#include "string.h"
#include "string"
#include <fstream>
#include <sstream>



HCPlot::HCPlot()
{
  cout << "<HCPlot::HCPlot>: LOADING" << endl;

  Es_max = (Variables -> gBeam) ;
  Ep_min = (Variables -> Ep_min);


  
}

HCPlot::~HCPlot()
{
  cout << "<HCPlot::~HCPlot>: END" << endl;
}

void HCPlot::SelectPlot()
{
  if(Variables -> plotPhaseSpace) Plot_phase_space();
  if(Variables -> plotStructureF) PlotSF();
  if( Variables -> plotCrossSecHe3) PlotHe3XS();

  
}



//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

void HCPlot:: Plot_phase_space()
{
  Es_min = (PSpace -> Es_min() );
  Ep_max = (PSpace -> Ep_max() );

  cout<<Es_max<<" "<<Es_min<<" "<<Ep_max<<" "<<Ep_min<<endl;

  //************** GRAPH THE PHASE-SPACE*****************

  TCanvas *c2;
  c2 = new TCanvas("c2","A Simple Graph Example",200,10,700,500);
  // c2->Divide(5,2);
  c2->cd();

  vgrEsDIS = PSpace -> vEs_PS;
  vgrEpDIS = PSpace -> vEp_PS;
  
  grDIS = new TGraph(vgrEsDIS.size(), &vgrEsDIS  [0] , &vgrEpDIS  [0]);
 
  // These values should be estimated from the kinematics
  // axisDIS = grDIS  ->GetXaxis();
  // axisDIS -> SetLimits(5. ,11. ); 
  //  grDIS   -> GetHistogram()->SetMaximum(8. );   // along          
  //  grDIS   -> GetHistogram()->SetMinimum(3. ); 
    
    
    
  grDIS -> SetFillColor(10);
  grDIS -> SetLineColor(kPink);
  grDIS -> SetLineWidth(4);
  grDIS -> SetMarkerColor(4);
  grDIS -> SetMarkerStyle(8);

  grDIS -> SetTitle(Form("Phase Space -- #Theta = %.2f deg", Variables->gTheta) );
  grDIS -> Draw("AL");
  grDIS -> GetHistogram()->SetXTitle("Es (GeV)");
  grDIS -> GetHistogram()->SetYTitle("Ep (GeV)");
      

  DISlineH = new TLine(Es_min,  Ep_min, Es_max, Ep_min  );
  DISlineH -> SetLineWidth(4);
  DISlineH -> SetLineColor(kPink);
  DISlineH -> Draw("same");
	   
	   
  DISlineV = new TLine(Es_max, Ep_min, Es_max, Ep_max   );
  DISlineV -> SetLineWidth(4);
  DISlineV -> SetLineColor(kPink);
  DISlineV -> Draw("same");
}


void HCPlot:: PlotSF()
{
  // The last Es bin, which is supposed to be Ebeam
  // from the construction, is the last vector in the F1/F2 vector
  lastbin  = Spectra -> NoEsBins;

  vnuPlot = Spectra -> vvnuCS[lastbin-1];
  //  vEpPlot = Spectra -> vvEpCS[lastbin-1];//the vector of Ep values
  // I don't think that copying the vector is efficient.
  // Maybe using <handler>->vector is more efficient

   
  if(Variables->nucleon == "proton")
    {
      nucle = "proton";  
      vvF1Plot = Model -> vvF1pro[lastbin-1];
      vvF2Plot = Model -> vvF2pro[lastbin-1];
      c3 = new TCanvas("c3","Structure Function",1);
      GraphSF();
    }
  else
    if(Variables->nucleon == "neutron")
    {
      nucle = "neutron";  
      vvF1Plot = Model -> vvF1neu[lastbin-1];
      vvF2Plot = Model -> vvF2neu[lastbin-1];
      c3 = new TCanvas("c3","Structure Function",1);
      GraphSF();
    }
  else
    if(Variables->nucleon == "He3")
      {
	c3 = new TCanvas("c3","Structure Function",1);
	c4 = new TCanvas("c4","Structure Function",1);
	c3->cd();
	nucle = "proton";  
	vvF1Plot = Model -> vvF1pro[lastbin-1];
	vvF2Plot = Model -> vvF2pro[lastbin-1];
	GraphSF();
	  
	c4->cd();
	nucle = "neutron";  
	vvF1Plot = Model -> vvF1neu[lastbin-1];
	vvF2Plot = Model -> vvF2neu[lastbin-1];
	GraphSF();

	
      }



  
}
void HCPlot::GraphSF()
{
  eBeam = Variables -> gBeam;

  // This function plots the last bin of the calculated Phase-Space previously
  // The last bin is supposed to be Es=Ebeam=10.38.

  // I think that if we declare the graph in the header, it will be cause some problems if we reuse the function
  TGraph *grF1;
  grF1 = new TGraph(vnuPlot.size(), &vnuPlot[0], &vvF1Plot[0]);	 
  
  grF1 -> SetTitle(Form("SF %s - %.2fdeg - %.2f GeV - W=%.2f GeV threshold", nucle, Variables-> gTheta, eBeam, Variables->WThreshold) );
  grF1 -> GetHistogram()->SetXTitle("#nu [GeV]");
  grF1 -> GetHistogram()->SetYTitle("F1/F2 ");

  
  if(Variables->plotStructureFrep=="line")
    {
      grF1 -> SetLineWidth(4);
      grF1 -> SetLineColor(kRed);
      grF1 -> Draw("AL");
    }
  else
    {
      grF1 -> SetMarkerColor(kBlue);
      grF1 -> SetMarkerStyle(7);
      grF1 -> Draw("AP");
    }
  
  // TCanvas *c4;
  //
  
  TGraph *grF2;
  grF2 = new TGraph(vnuPlot.size(), &vnuPlot[0], &vvF2Plot[0]);	 
  
  if(Variables->plotStructureFrep=="line")
    {
      grF2 -> SetLineWidth(4);
      grF2 -> SetLineColor(kBlue);
      grF2 -> Draw("SAMEL");
    }
  else
    {
      grF2 -> SetMarkerColor(kBlue);
      grF2 -> SetMarkerStyle(7);
      grF2 -> Draw("SAMEP");
    }
  
   auto legend = new TLegend(0.15,0.7,0.3,0.9);
   legend->SetFillStyle(0);
   legend->SetBorderSize(0);
   legend->AddEntry(grF1,"F1","l");
   legend->AddEntry(grF2,"F2","l");
   legend->Draw();
    

  return;
}


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
void HCPlot::PlotHe3XS()
{
  lastbin = Spectra -> NoEsBins;
  vnuPlot = Spectra -> vvnuCS[lastbin-1];
 
  vXS = Model -> vvXS[lastbin-1];
  
  cHe3XS = new TCanvas("cHeXS","He3 Cross Section",1);

  grXS = new TGraph(vnuPlot.size(), &vnuPlot[0], &vXS[0]);	 

      
  TAxis *axisINE; 
  grXS  ->SetLineWidth(4);
  grXS ->SetMarkerColor(4);
  grXS ->SetLineColor(kRed);
  
  // axisINE =grXS ->GetXaxis();
  // axisINE -> SetLimits(20.,180.); 
  //grXS -> GetHistogram()->SetMinimum(0.);

  grXS -> SetTitle(Form("KIN A - %.2fdeg - 10.38GeV - W=%.2f GeV threshold", Variables->gTheta, Variables->WThreshold ));
  grXS -> GetHistogram()->SetXTitle("#nu [GeV]");
  grXS -> GetHistogram()->SetYTitle("Cross-Section [nb/MeV/sr]");
  grXS -> Draw("AL");

  

}

HCPlot *Plot = NULL;
