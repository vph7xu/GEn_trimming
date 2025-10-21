
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is to take analyzed data and
//   to calculate the asymmetry. It requires the output root
//   file from QuasiElastic_ana.C.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../include/gen-ana.h"


void Asymmetry_yield(const char *configfilename,std::string filebase="../outfiles/QE_data")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // seting up the desired SBS configuration
  TString conf = jmgr->GetValueFromKey_str("GEN_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromKey_str("Ntype");

  int model = jmgr->GetValueFromKey<int>("model");

  int IHWP_Flip = jmgr->GetValueFromKey<int>("IHWP_Flip");

  // elastic cut limits
  double W2min = jmgr->GetValueFromKey<double>("W2min");
  double W2max = jmgr->GetValueFromKey<double>("W2max");

  double dy_min = jmgr->GetValueFromKey<double>("dymin");
  double dy_max = jmgr->GetValueFromKey<double>("dymax");


  TString inFile = Form("%s_" + sbsconf.GetSBSconf() + "_sbs%dp_nucleon_%s_model%d.root", 
			 filebase.c_str(),   sbsconf.GetSBSmag(), Ntype.c_str(), model);
  TFile *fin = new TFile(inFile.Data(), "read");
  TTree *T = (TTree*)fin->Get("Tout");

  T->SetBranchStatus("*",0);

  // This is the variables we need from the analyzed root file
  int runnum;   setrootvar::setbranch(T,"runnum","",&runnum);
  bool WCut;   setrootvar::setbranch(T,"WCut","",&WCut);
  bool pCut;   setrootvar::setbranch(T,"pCut","",&pCut);
  bool nCut;   setrootvar::setbranch(T,"nCut","",&nCut);
  bool coinCut;   setrootvar::setbranch(T,"coinCut","",&coinCut);
  double W2;   setrootvar::setbranch(T,"W2","",&W2);
  double dx;   setrootvar::setbranch(T,"dx","",&dx);
  double dy;   setrootvar::setbranch(T,"dy","",&dy);
  double coin_time;   setrootvar::setbranch(T,"coinT_trig","",&coin_time);
  double hcal_time;   setrootvar::setbranch(T,"hcal_time","",&hcal_time);
  double hodo_time[1000];   setrootvar::setbranch(T,"hodo_time","",&hodo_time);
  int helicity;   setrootvar::setbranch(T,"helicity","",&helicity);
  int IHWP;   setrootvar::setbranch(T,"IHWP","",&IHWP);

  
  TH1F *hdx_p = new TH1F("hdx_p","#Deltax for helicity +1;#Deltax;Entries",100,-6,4);
  TH1F *hdx_m = new TH1F("hdx_m","#Deltax for helicity -1;#Deltax;Entries",100,-6,4);

  //Totals for all runs
  int Yp_n_total = 0;
  int Ym_n_total = 0;
  int ncut_n_total = 0;
  int Yp_p_total = 0;
  int Ym_p_total = 0;
  int ncut_p_total = 0;

  //Y+ and Y- neutron helicity yields
  //For individual runs
  int Yp_n = 0;
  int Ym_n = 0;
  int ncut_n = 0;
  
  //Y+ and Y- proton helicity yields
  //For individual runs
  int Yp_p = 0;
  int Ym_p = 0;
  int ncut_p = 0;

  // Setup some variables for helicity per run
  vector<double> runs;
  vector<double> A_n;
  vector<double> A_n_err;
  vector<double> A_p;
  vector<double> A_p_err;
  vector<double> p_nevents_points, n_nevents_points, A_n_nevents, A_p_nevents, A_n_err_nevents, A_p_err_nevents;

  int nevent = 0;
  int nQE = 0;

  //This will loop over all events and record the helicity
  while(T->GetEntry(nevent++)){
    
    if(IHWP == 1) helicity *= -1*IHWP_Flip;
    else if(IHWP == -1) helicity *= 1*IHWP_Flip;
    else continue; //skip if the helicity is undefined
    
    //Cuts for good elastic events
    if(W2 > W2min && W2 < W2max && coinCut && dy > dy_min && dy < dy_max){
      nQE++;
      if(helicity == 1)
	hdx_p->Fill(dx);
      if(helicity == -1)
	hdx_m->Fill(dx);

    }
  }

  //Plot dx for positive helicity
  TCanvas *c1 = new TCanvas("c1","",800,600);
  hdx_p->Draw();
  TF1 *fit_n;
  TF1 *fit_p;
  TF1 *fit_bg;
  Analysis::He3_fit(hdx_p,conf,&fit_bg,&fit_n,&fit_p);
  fit_n->Draw("same");
  //fit_p->Draw("same");
  fit_bg->SetLineStyle(2);
  fit_bg->Draw("same");

  double n_p_yield = abs(fit_n->Integral(-4,4)/hdx_p->GetBinWidth(0));
  double p_p_yield = abs(fit_p->Integral(-5,0)/hdx_p->GetBinWidth(0));

  //Plot dx for negative helicity
  TCanvas *c2 = new TCanvas("c2","",800,600);
  hdx_m->Draw();
  Analysis::He3_fit(hdx_m,conf,&fit_bg,&fit_n,&fit_p);
  fit_n->Draw("same");
  //fit_p->Draw("same");
  fit_bg->SetLineStyle(2);
  fit_bg->Draw("same");

  double n_m_yield = abs(fit_n->Integral(-4,4)/hdx_m->GetBinWidth(0));
  double p_m_yield = abs(fit_p->Integral(-5,0)/hdx_m->GetBinWidth(0));
  double n_total = (n_p_yield + n_m_yield);
  double p_total = (p_p_yield + p_m_yield);

  double A_n_total = (n_p_yield - n_m_yield)*1.0/n_total;
  double Err_n = sqrt((1 - A_n_total*A_n_total)/n_total);

  double A_p_total = (p_p_yield - p_m_yield)*1.0/p_total;
  double Err_p = sqrt((1 - A_p_total*A_p_total)/p_total);

  cout<<"Total Neutron Events "<<n_total<<endl;
  cout<<"Y+ "<<n_p_yield<<endl;
  cout<<"Y- "<<n_m_yield<<endl;
  cout<<"A = "<<A_n_total*100<<"% +/- "<<Err_n*100<<"%"<<endl;
  cout<<endl;
  cout<<"Total Proton Events "<<p_total<<endl;
  cout<<"Y+ "<<p_p_yield<<endl;
  cout<<"Y- "<<p_m_yield<<endl;
  cout<<"A = "<<A_p_total*100<<"% +/- "<<Err_p*100<<"%"<<endl;
}
