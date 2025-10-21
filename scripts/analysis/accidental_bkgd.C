//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified February 25, 2024
//
//
//   The purpose of this script is to use He3 data to test
//   cuts for accidental background contributions.
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


void accidental_bkgd(TString cfg){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TString DB_file = "../../DB/Helicity_quality.csv";
  getDB(DB_file);

  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file);

  // Set variables for cuts
  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;

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

  double coin_bg_low = kin_info.coin_time_cut[0] + (1 + 3)*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  double coin_bg_hi = kin_info.coin_time_cut[0] + (3 + 3)*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];

  // Set cuts limits for a number of points
  const int N_points = 20;
  double coin_low[N_points], coin_hi[N_points];
  int N_p[N_points], N_m[N_points];

  for(int i=0; i < N_points; i++){
    coin_low[i] = coin_bg_low + i*kin_info.coin_time_cut[1];
    coin_hi[i] = coin_low[i] + 2*kin_info.coin_time_cut[1];
    N_p[i] = 0;
    N_m[i] = 0;
  }

  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_info,1,1);

  TH1F *hcoin_time = new TH1F("hcoin_time","Coincidence Time;Coincidence Time (ns);Entries",100,40,180);

  int nevent = 0;
  int maxevent = T_data->fChain->GetEntries();
  int total = 0;

  while(nevent < maxevent){
    T_data->GetEntry(nevent++);   

    // Use all QE cuts 
    if(T_data->helicity != -1 && T_data->helicity != 1) continue;
    if(T_data->W2 < W2min || T_data->W2 > W2max) continue;
    if(!fGoodHel[T_data->runnum]) continue;
    if(T_data->dy < dymin || T_data->dy > dymax) continue;
    if(T_data->dx < dxmin || T_data->dx > dxmax) continue;

    int helicity = T_data->helicity;
    helicity *= -1*T_data->IHWP*IHWP_Flip;
    
    hcoin_time->Fill(T_data->coin_time);  // Fill coincidence histogram
    
    for(int i=0; i < N_points; i++){
      if(T_data->coin_time > coin_low[i] && T_data->coin_time < coin_hi[i]){
	if(helicity == 1){
	  N_p[i]++;
	}
	else if(helicity == -1){
	  N_m[i]++;
	}
      }
    }

  }

  // print helicity for each cut test
  for(int i=0; i < N_points; i++){
    cout<<i<<" "<<N_p[i]<<" "<<N_m[i]<<" "<<1.0*(N_p[i] - N_m[i]) / (N_p[i] + N_m[i])<<endl;
  }

  TCanvas *c = new TCanvas("c","",800,600);
  hcoin_time->Draw();

  TLine *lline = new TLine(coin_bg_low,0,coin_bg_low,hcoin_time->GetMaximum());
  lline->SetLineColor(kRed);
  lline->Draw("same");

  TLine *rline = new TLine(coin_bg_hi,0,coin_bg_hi,hcoin_time->GetMaximum());
  rline->SetLineColor(kRed);
  rline->Draw("same");

  TPaveText *pt = new TPaveText(.65,.50,.88,.70,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText(Form("%g < W^{2} < %g",kin_info.W2min,kin_info.W2max));
  pt->AddText(Form("|#DeltaR| < %g",dxmax));
  pt->SetFillColor(0);
  pt->Draw("same");
 
}


