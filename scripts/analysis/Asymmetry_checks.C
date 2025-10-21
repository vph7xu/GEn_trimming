//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified February 2, 2024
//
//
//   The purpose of this script is to compare real data and
//   simulated data for the same kinematic point
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

void Asymmetry_checks(TString cfg = "GEN2",TString tgt = "He3"){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TString DB_file = "../../DB/Helicity_quality.csv";
  getDB(DB_file);

  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file);

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

  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_info,1,1);
  
  /////Set the histograms
  int nbinsdx = 100;
  double xmin = -4;
  double xmax = 2.5;
  double W2min_r = -1;
  double W2max_r = 5;
  double W2bin_size = 0.2;

  if(cfg == "GEN2"){
    xmin = -6;
    xmax = 3;
    W2min_r = -1;
    W2max_r = 5;
  }
  else if(cfg == "GEN3"){
    xmin = -6;
    xmax = 3;
    W2min_r = -1;
    W2max_r = 6;
  }
  else if(cfg == "GEN4"){
    xmin = -6;
    xmax = 3;
    W2min_r = -2;
    W2max_r = 6;
  }

  const int nbinsW2 = (int)(W2max_r - W2min_r) / W2bin_size;
  
  
  //dx
  TH2F *hW2_x = new TH2F("hW2_x","",nbinsW2,W2min_r,W2max_r,nbinsdx,xmin,xmax);
  
  double W2_array[nbinsW2];
  double A_array[nbinsW2];
  double A_err_array[nbinsW2];
  int Yp_array[nbinsW2];
  int Ym_array[nbinsW2];

  for(int ibin=0; ibin < nbinsW2; ibin++){
    Yp_array[ibin] = 0;
    Ym_array[ibin] = 0;
  }
  
  int nevent = 0;
  int maxevent = T_data->fChain->GetEntries();
  int N_p = 0, N_m = 0;
  
  while(nevent < maxevent){
    T_data->GetEntry(nevent++);
    
    int helicity = T_data->helicity;
    helicity *= -1*T_data->IHWP*IHWP_Flip; 

    if(!fGoodHel[T_data->runnum]) continue;   
    if(helicity != -1 && helicity != 1) continue;
    if(T_data->coin_time < coin_min || T_data->coin_time > coin_max) continue;
    if(T_data->W2 < W2min_r || T_data->W2 > W2max_r) continue;
    if(T_data->dy > dy_bg_min && T_data->dy < dy_bg_max) continue;

    hW2_x->Fill(T_data->W2,T_data->dx);
 
    int W2bin = (int) ((T_data->W2 - W2min_r) / W2bin_size);

    if(T_data->dx > dxmin && T_data->dx < dxmax){ //Cut around neutron in dx
      if(helicity == 1){
	Yp_array[W2bin]++;
	N_p++;
      }
      else if(helicity == -1){
	Ym_array[W2bin]++;
	N_m++;
      }
    }
  }
  
  for(int ibin=0; ibin < nbinsW2; ibin++){
    W2_array[ibin] = W2min_r + (ibin + 1) * W2bin_size - W2bin_size/2;
    A_array[ibin] = (Yp_array[ibin] - Ym_array[ibin])*1.0/(Yp_array[ibin] + Ym_array[ibin]);
    A_err_array[ibin] = sqrt((1 - A_array[ibin]*A_array[ibin])/(Yp_array[ibin] + Ym_array[ibin])) * 100;
    A_array[ibin] *= 100; //convert to percentage
  }

  TGraphErrors *gA = new TGraphErrors(nbinsW2,W2_array,A_array,0,A_err_array);
 
  gA->SetTitle(";W^{2} (GeV^{2});A (%)");
  hW2_x->SetTitle(cfg + " Kinematics with Asymmetry;;#Deltax (m)");
  gA->SetMarkerStyle(8);

  gA->GetXaxis()->SetLabelSize(0.1);
  gA->GetYaxis()->SetLabelSize(0.1);
  hW2_x->GetYaxis()->SetLabelSize(0.05);
  gA->GetXaxis()->SetTitleSize(0.10);
  gA->GetYaxis()->SetTitleSize(0.10);
  hW2_x->GetYaxis()->SetTitleSize(0.05);
  gA->GetYaxis()->SetTitleOffset(0.28);
  hW2_x->GetYaxis()->SetTitleOffset(0.35);

  TCanvas *c = new TCanvas("c","",800,600);
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  
  pad1->SetBottomMargin(0.00); // Set bottom margin for pad1
  pad2->SetTopMargin(0.00); // Set top margin for pad2
  pad2->SetBottomMargin(0.30);
  
  pad1->Draw();
  pad2->Draw();

  TLine *lymin = new TLine(W2min_r,dxmin,W2max_r,dxmin);
  TLine *lymax = new TLine(W2min_r,dxmax,W2max_r,dxmax);
  lymin->SetLineColor(kRed);
  lymax->SetLineColor(kRed);
  
  pad1->cd();
  pad1->SetGridx();
  hW2_x->Draw("colz");
  lymin->Draw("same");
  lymax->Draw("same");

  TPaveText *pt = new TPaveText(.13,.10,.38,.20,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText("Coincidence Cuts");
  pt->SetFillColor(0);
  pt->Draw("same");

  TLegend *legend = new TLegend(0.13,0.03,0.38,0.10);
  legend->AddEntry(lymin,"#Deltax cut applied below","l");
  legend->SetLineColor(0);
  legend->Draw("same");
  
  pad2->cd();
  pad2->SetGridx();
  gA->Draw("AP");
  gA->GetXaxis()->SetLimits(W2min_r,W2max_r);
  gA->GetYaxis()->SetRangeUser(-3,8);
  
  TLine *l0 = new TLine(W2min_r,0,W2max_r,0);
  l0->Draw("same");

  cout<<N_p - N_m<<" "<<N_p + N_m<<" "<<1.0*(N_p - N_m) / (N_p + N_m)<<endl;

  TString output = "Asymmetry_W2bins_"+cfg+".pdf";
  
  c->SaveAs("../../plots/" + output);
  
}
