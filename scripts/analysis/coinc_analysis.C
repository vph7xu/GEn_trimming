//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is to calculate the coincidence
//   time between the BB and SBS arms. It will then print out
//   timing cuts that should be used for coincidence.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"


void do_fit_pol4(TCanvas *c, TH1D *h1,double bg_low, double bg_high, double signal_low,double signal_high){

  c->cd();
  double par[8];
  TF1 *gaus_func = new TF1("gaus_func","gaus",signal_low,signal_high);
  TF1 *bg_func = new TF1("bg_func","pol4",bg_low,bg_high);
  TF1 *total_func = new TF1("total_func","gaus(0) + pol4(3)",bg_low,bg_high);

  h1->Fit("gaus_func","qNR");
  h1->Fit("bg_func","qNR+");

  gaus_func->GetParameters(&par[0]);
  bg_func->GetParameters(&par[3]);

  total_func->SetParameters(par);
  h1->Fit(total_func,"qNR+"); 
  total_func->GetParameters(&par[0]);
  //gaus_func->SetParameters(&par[0]);

  total_func->Draw("same");

  //cout<<par[1]<<" +/- "<<par[2]<<endl;

  TPaveText *pt = new TPaveText(.11,.6,.45,.73,"ndc");
  pt->AddText(Form("Coin Time = %g +/- %g ns",par[1],par[2]));
  pt->SetFillColor(0);
  pt->Draw("same");

}

void do_fit_pol1(TCanvas *c, TH1D *h1,double bg_low, double bg_high, double signal_low,double signal_high){

  c->cd();
  double par[5];
  TF1 *gaus_func = new TF1("gaus_func","gaus",signal_low,signal_high);
  TF1 *bg_func = new TF1("bg_func","pol1",bg_low,bg_high);
  TF1 *total_func = new TF1("total_func","gaus(0) + pol1(3)",bg_low,bg_high);

  h1->Fit("gaus_func","qNR");
  h1->Fit("bg_func","qNR+");
  
  gaus_func->GetParameters(&par[0]);
  bg_func->GetParameters(&par[3]);
  
  total_func->FixParameter(3,par[3]);
  total_func->FixParameter(4,par[4]);
  total_func->SetParLimits(0,par[0]*0.95,par[0]*1.05);
  total_func->SetParLimits(1,par[1]*0.95,par[1]*1.05);
  total_func->SetParLimits(2,par[2]*0.95,par[2]*1.05);

  total_func->SetParameters(par);
  h1->Fit(total_func,"qNR+"); 
  //total_func->GetParameters(&par[0]);

  total_func->Draw("same");
   
  TPaveText *pt = new TPaveText(.11,.6,.45,.73,"ndc");
  pt->AddText(Form("Coin Time = %g +/- %g ns",par[1],par[2]));
  pt->SetFillColor(0);
  pt->Draw("same");

}


void coinc_analysis(TString cfg = "GEN2"){

  gStyle->SetOptStat(0);

  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";
  Utilities::KinConf He3_kin = Utilities::LoadKinConfig(jmgr_file);

  jmgr_file = "../../config/" + cfg + "_H2.cfg";
  Utilities::KinConf H2_kin = Utilities::LoadKinConfig(jmgr_file);

  analyzed_tree *T_He3 = Utilities::LoadAnalyzedRootFiles(He3_kin,1,1);
  analyzed_tree *T_H2 = Utilities::LoadAnalyzedRootFiles(H2_kin,0,1);

  vector<double> dx_n = He3_kin.dx_n;
  double Nsigma_dx_n = He3_kin.Nsigma_dx_n;
  vector<double> dy_n = He3_kin.dy_n;
  double Nsigma_dy_n = He3_kin.Nsigma_dy_n;
  double dxmin = dx_n[0] - dx_n[1];
  double dxmax = dx_n[0] + dx_n[1];
  double dymin = dy_n[0] - dy_n[1];
  double dymax = dy_n[0] + dy_n[1];

  TH1D *h_coinc_H2_cal = new TH1D("h_coinc_H2_cal","H2 Coincidence Time;HCal TDC Time (ns) - BBCal ADC Time; Entries",150,-10,200);
  TH1D *h_coinc_H2_hodo = new TH1D("h_coinc_H2_hodo","H2 Coincidence Time;HCal TDC Time (ns) - Hodo TDC Time; Entries",150,-10,200);
  
  int nevent = 0;
  int maxevent = T_H2->fChain->GetEntries();
  
  //Loop over all events on the H2 file
  while(nevent < maxevent){
    T_H2->GetEntry(nevent++); 

    if(T_H2->W2 < H2_kin.W2min || T_H2->W2 > H2_kin.W2max) continue;
    if(T_H2->dy < H2_kin.dymin || T_H2->dy > H2_kin.dymax) continue;

    h_coinc_H2_cal->Fill(T_H2->hcal_time - T_H2->bbcal_time);
    h_coinc_H2_hodo->Fill(T_H2->hcal_time - T_H2->hodo_time[0]);
    
  }
  
  
  TH1D *h_coinc_He3_cal = new TH1D("h_coinc_He3_cal","He3 Coincidence Time;HCal ADC Time (ns) - BBCal ADC Time; Entries",150,40,180);
  TH1D *h_coinc_He3_hodo = new TH1D("h_coinc_He3_hodo","He3 Coincidence Time;HCal TDC Time (ns) - Hodo TDC Time; Entries",150,-100,100);
  

  nevent = 0;
  maxevent = T_He3->fChain->GetEntries();

  //Loop over all events on the H2 file
  while(nevent < maxevent){
    T_He3->GetEntry(nevent++); 

    if(T_He3->W2 < He3_kin.W2min || T_He3->W2 > He3_kin.W2max) continue;
    if(T_He3->dy < dymin || T_He3->dy > dymax) continue;
    if(T_He3->dx < dxmin || T_He3->dx > dxmax) continue;
    
    h_coinc_He3_cal->Fill(T_He3->hcal_time - T_He3->bbcal_time);
    h_coinc_He3_hodo->Fill(T_He3->hcal_time - T_He3->hodo_time[0]); 
    
  }

  TCanvas *c1 = new TCanvas("c1","",800,600);
  h_coinc_H2_cal->Draw();
  do_fit_pol4(c1,h_coinc_H2_cal,-10,200,60,100); //Fit H2 data

  TCanvas *c2 = new TCanvas("c2","",800,600);
  h_coinc_He3_cal->Draw();
  do_fit_pol1(c2,h_coinc_He3_cal,40,200,89,95); //Fit He3 data


}
