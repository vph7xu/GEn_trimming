//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified February 29, 2024
//
//
//   The purpose of this script is to calculate the helicity
//   of inelastic data.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"

map<int,double> fGoodHel;


void getDB(TString DB_file){

  // Entries should follow this form:
  //{var name, variable pointer, vairable description, 1/0 (mandatory/not mandatory variable)}
  DBparse::DBRequest request[] = {
    {"Good Helicity", &fGoodHel, NULL, "Is the helicity readback good (0/1 = no/yes)", 1}
  };
  
  const int nvar = sizeof(request) / sizeof(request[0]);
  
  DB_load(DB_file,request,nvar);
}


void Inelastic_contamination(TString cfg){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TString DB_file = "../../DB/Helicity_quality.csv";
  getDB(DB_file);

  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file);

   // elastic cut limits
  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;

  double dy_bg_min = kin_info.dymin;
  double dy_bg_max = kin_info.dymax;

  vector<double> dx_n = kin_info.dx_n;
  double Nsigma_dx_n = kin_info.Nsigma_dx_n;
  vector<double> dy_n = kin_info.dy_n;
  double Nsigma_dy_n = kin_info.Nsigma_dy_n;
  double dxmin = dx_n[0] - dx_n[1];
  double dxmax = dx_n[0] + dx_n[1];
  double dymin = dy_n[0] - dy_n[1];
  double dymax = dy_n[0] + dy_n[1];

  int IHWP_Flip = kin_info.IHWP_Flip;

  double coin_min = kin_info.coin_time_cut[0] - kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  double coin_max = kin_info.coin_time_cut[0] + kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  
  double coin_bg_min = kin_info.coin_time_cut[0] + (1 + 3)*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  double coin_bg_max = kin_info.coin_time_cut[0] + (3 + 3)*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];

  SetHe3Pol();

  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_info,1,0);
  analyzed_tree *T_sim = Utilities::LoadAnalyzedRootFiles(kin_info,0,0);

  distribution_fits *dists = new distribution_fits();

  //if(cfg == "GEN2") dists->SetBgShapeOption("pol4");
  //else dists->SetBgShapeOption("from data");
  dists->SetBgShapeOption("from data");


  /////Set the histograms
  int nbins = 100;
  double hxmin = -4;
  double hxmax = 2.5;

  if(cfg == "GEN2"){
    hxmin = -6;
    hxmax = 3;
  }
  
  //dx
  TH1F *hdx_data = new TH1F("hdx_data","",nbins,hxmin,hxmax);
  TH1F *hdx_sim_p = new TH1F("hdx_sim_p","",nbins,hxmin,hxmax);
  TH1F *hdx_sim_n = new TH1F("hdx_sim_n","",nbins,hxmin,hxmax);
  TH1F *hdx_bg_data = new TH1F("hdx_bg_data","",nbins,hxmin,hxmax);

  //W2
  TH1F *hW2_all = new TH1F("hW2_all","W^{2} With Different Cuts;W^{2};",nbins,-1,7); // All W2 in the dx cuts
  TH1F *hW2_cut1 = new TH1F("hW2_cut1","",nbins,-1,7); // W2 in inelastic cuts
  TH1F *hW2_cut2 = new TH1F("hW2_cut2","",nbins,-1,7); // W2 in elastic cuts
  
  TCut CutSimP = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 1) * weight",W2min,W2max,dymin,dymax);
  TCut CutSimN = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 0) * weight",W2min,W2max,dymin,dymax);
  
  T_sim->fChain->Draw("dx>>hdx_sim_p",CutSimP);
  T_sim->fChain->Draw("dx>>hdx_sim_n",CutSimN);

  int nevent = 0;
  int maxevent = T_data->fChain->GetEntries();
  int total = 0;
  int N_p =0, N_m = 0;
  analyzed_info *Asym_total = new analyzed_info();

  while(nevent < maxevent){
    T_data->GetEntry(nevent++);  

    int helicity = T_data->helicity;
    helicity *= -1*T_data->IHWP*IHWP_Flip; 

    if(!fGoodHel[T_data->runnum]) continue;
    if(T_data->helicity != -1 && T_data->helicity != 1) continue;
    if(T_data->coin_time < coin_min || T_data->coin_time > coin_max) continue;

    // Here we will different histograms for different cuts
    if(T_data->dx > dxmin && T_data->dx < dxmax){
      hW2_all->Fill(T_data->W2);
      
      // We count inelastics outside a wide dy region
      if(T_data->dy < dy_bg_min || T_data->dy > dy_bg_max){	
	hW2_cut1->Fill(T_data->W2);
	Asym_total->IterateInelasticCount(helicity); // Add helicity counts
      }
      if(T_data->dy > dxmin && T_data->dy < dxmax){
	hW2_cut2->Fill(T_data->W2);
      }
    }

    if(T_data->W2 < W2min || T_data->W2 > W2max) continue;
    if(T_data->dy < dy_bg_min || T_data->dy > dy_bg_max) hdx_bg_data->Fill(T_data->dx);
    if(T_data->dy < dymin || T_data->dy > dymax) continue;

    // Fill elastic histograms as usual
    hdx_data->Fill(T_data->dx);

    if(T_data->dx > dxmin && T_data->dx < dxmax){ //Cut around neutron in dx

      UpdateExpansionCoefficients(T_data);      

      if(helicity == 1){
	N_p++;
      }
      else if(helicity == -1){
	N_m++;
      }
    }
  }
  
  
  hW2_cut1->SetLineColor(kRed);
  hW2_cut2->SetLineColor(kGreen);

  TCanvas *c = new TCanvas("c","",800,600);
  hW2_all->Draw();
  hW2_cut1->Draw("same");
  hW2_cut2->Draw("same");

  TLegend *legend = new TLegend(0.11,0.72,0.45,0.89);
  legend->AddEntry("hW2_all",Form("|#Deltax| < %g",dxmax),"l");
  legend->AddEntry("hW2_cut1",Form("|#Deltax| < %g & |#Deltay| > %g",dxmax,dy_bg_max),"l");
  legend->AddEntry("hW2_cut2",Form("|#DeltaR| < %g",dxmax),"l");
  legend->SetLineColor(0);
  legend->Draw("same");
 
  dists->SetDataShape(hdx_data);
  dists->SetPShape(hdx_sim_p);
  dists->SetNShape(hdx_sim_n);
  dists->SetBgShape(hdx_bg_data);
  
  dists->He3_fit_dists();

  int in_yield = dists->GetBgYield(dxmin, dxmax);

  double Delta_in = Asym_total->N_in_p - Asym_total->N_in_m;
  double Sigma_in = Asym_total->N_in_p + Asym_total->N_in_m;
  double Sigma_inel = in_yield;
  double Sigma = N_p + N_m;
  double A_in_old = Delta_in / Sigma_in;
  double A_in = Delta_in * Sigma_inel / (Sigma * Sigma_in);

  cout<<Asym_total->N_in_p<<" "<<Asym_total->N_in_m<<" "<<in_yield<<endl;
  cout<<A_in_old<<" "<<A_in<<endl;

  //Copy all the result histograms
  TH1F *hdx_data_plot = dists->GetDataHist();
  TH1F *hdx_sim_p_plot = dists->GetPHist();
  TH1F *hdx_sim_n_plot = dists->GetNHist();
  TH1F *hdx_bg_plot = dists->GetBgHist();
  TH1F *hdx_total_fit_plot = dists->GetTotalHist();

  hdx_data->SetMarkerStyle(kFullCircle);
  hdx_total_fit_plot->SetFillColorAlpha(30,0.5);
  hdx_sim_p->SetFillColorAlpha(kRed,0.3);
  hdx_sim_n->SetFillColorAlpha(kBlue,0.3);
  hdx_bg_plot->SetFillColorAlpha(kMagenta,0.3);

  hdx_total_fit_plot->SetLineStyle(7);
  hdx_sim_p->SetLineStyle(7);
  hdx_sim_n->SetLineStyle(7);
  hdx_bg_plot->SetLineStyle(7);
  
  hdx_total_fit_plot->SetLineColor(30);
  hdx_sim_p->SetLineColor(kRed);
  hdx_sim_n->SetLineColor(kBlue);
  hdx_bg_plot->SetLineColor(kMagenta);
  
  TCanvas *c2 = new TCanvas("c2","",800,600);
  hdx_data->Draw();
  hdx_data_plot->Draw("hist");
  hdx_sim_p_plot->Draw("same hist");
  hdx_sim_n_plot->Draw("same hist");
  hdx_bg_plot->Draw("same hist");

}


