//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified March 4, 2024
//
//
//   The purpose of this script is to calculate the pion 
//   asymmetry contribution
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"

TH1F *hpse_sim_He3;
TH1F *hpse_sim_pim;

// Fit for preshower simulation
double fitdist( double *x, double *par){
  double dx = x[0];
  
  double Norm_He3 = par[0];
  double Norm_pi = par[1];
  double ps_scale_He3 = par[2];
  double ps_scale_pi = par[3];
    
  double simu = Norm_He3 * hpse_sim_He3->Interpolate(ps_scale_He3 * dx) + Norm_pi * hpse_sim_pim->Interpolate(ps_scale_pi * dx);
  
  return simu;   
}

void pion_contamination(){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  TString jmgr_He3 = "../../config/GEN2_He3.cfg";
  TString jmgr_pim = "../../config/GEN2_pim_sim.cfg";
  Utilities::KinConf kin_He3 = Utilities::LoadKinConfig(jmgr_He3);
  Utilities::KinConf kin_pim = Utilities::LoadKinConfig(jmgr_pim);
  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_He3,1,1);
  analyzed_tree *T_sim_He3 = Utilities::LoadAnalyzedRootFiles(kin_He3,0,1);
  analyzed_tree *T_sim_pim = Utilities::LoadAnalyzedRootFiles(kin_pim,0,1);

  
  int nbins = 150;
  double hxmin = 0;
  double hxmax = 2;

  TH1F *hpse_data = new TH1F("hpse_data","Pion Data vs Simulation;Preshower Energy (GeV);",nbins,hxmin,hxmax);
  hpse_sim_He3 = new TH1F("hpse_sim_He3","",nbins,hxmin,hxmax);
  hpse_sim_pim = new TH1F("hpse_sim_pim","",nbins,hxmin,hxmax);
  TH1F *hpse_sim_tot = new TH1F("hpse_sim_tot","",nbins,hxmin,hxmax);
  
  TCut simNCut = "(ePS > 0.01 && WCut && fnucl == 0)*weight";
  TCut simPiCut = "(ePS > 0.01 && WCut)*weight";
  
  T_data->fChain->Draw("ePS>>hpse_data","coinCut && WCut && nCut");
  T_sim_He3->fChain->Draw("ePS>>hpse_sim_He3",simNCut);
  T_sim_pim->fChain->Draw("ePS>>hpse_sim_pim",simPiCut);

  double scale_data = 1.0/hpse_data->Integral();
  double scale_He3 = 1.0/hpse_sim_He3->Integral();
  double scale_pim = 1.0/hpse_sim_pim->Integral();

  hpse_data->Scale(scale_data);
  hpse_sim_He3->Scale(scale_He3);
  hpse_sim_pim->Scale(scale_pim);

  //Set fit to function fitsim
  TF1 *FitFunc = new TF1( "FitFunc", fitdist,hxmin,hxmax,4);

  //Set some arbitrary starting values, should not be hardcoding this
  FitFunc->SetNpx(1000);
  double startpar[] = {1.0,1.0,1.0,1.0};
  FitFunc->SetParameters(startpar);
  FitFunc->SetParLimits(0,0.01,1);
  FitFunc->SetParLimits(1,0.01,1);
  FitFunc->SetParLimits(2,0.01,1);
  FitFunc->SetParLimits(3,0.01,1);

  hpse_data->Fit(FitFunc,"0","",hxmin,hxmax);

  FitFunc->SetParameter(0,0.92);
  FitFunc->SetParameter(1,0);
  FitFunc->SetParameter(2,0.86);
  FitFunc->SetParameter(3,1);

  hpse_sim_He3->Reset();
  hpse_sim_pim->Reset();

  T_sim_He3->fChain->Draw(Form("%g*ePS>>hpse_sim_He3",FitFunc->GetParameter(2)),simNCut);
  T_sim_pim->fChain->Draw(Form("%g*ePS>>hpse_sim_pim",FitFunc->GetParameter(3)),simPiCut);
  
  hpse_data->Scale(1.0/scale_data);
  hpse_sim_He3->Scale(scale_He3*FitFunc->GetParameter(0)/scale_data);
  hpse_sim_pim->Scale(scale_pim*FitFunc->GetParameter(1)/scale_data);

  for(int ibin = 0; ibin < nbins;ibin++){
    hpse_sim_tot->SetBinContent(ibin,hpse_sim_He3->GetBinContent(ibin) + hpse_sim_pim->GetBinContent(ibin));
  }

  hpse_sim_He3->SetLineColor(kRed);
  hpse_sim_pim->SetLineColor(kGreen);
  hpse_sim_tot->SetLineColor(kMagenta);

  TCanvas *c = new TCanvas("c","",800,600);
  hpse_data->Draw("hist");
  hpse_sim_pim->Draw("same hist");    
  hpse_sim_He3->Draw("same hist");
  //hpse_sim_tot->Draw("same hist");
  
  TLegend *legend = new TLegend(0.65,0.72,0.89,0.89);
  legend->AddEntry("hpse_data","He3 Data","l");
  legend->AddEntry("hpse_sim_He3","QE He3 Sim","l");
  legend->AddEntry("hpse_sim_pim","#pi^{-} Sim","l");
  legend->AddEntry("hpse_sim_tot","Sim Total","l");
  legend->SetLineColor(0);
  legend->Draw("same");
 
}


