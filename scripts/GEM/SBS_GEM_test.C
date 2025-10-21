
#include "../../include/gen-ana.h"

void SBS_GEM_test(){

  double hcal_dist = 12; //Dist from hcal to GEM
  double z0_start = 5.0;
  double offset_start = 0.045;

  //TString rootdir = "/lustre19/expphy/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/Rootfiles/GEN4/He3/rootfiles";
  TString rootdir = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/TEST/rootfiles";

  TChain *T = new TChain("T");
  T->Add(rootdir + "/*6083*");

  T->SetBranchStatus("*",0);

  int maxhits = 1000;

  double bb_ntrack, sbs_ntrack;
  double bb_vz[maxhits], bb_tgy[maxhits], bb_px[maxhits], bb_py[maxhits], bb_pz[maxhits], bb_p[maxhits], bb_th[maxhits];
  double sbs_ph[maxhits], sbs_th[maxhits], sbs_y[maxhits], sbs_x[maxhits];
  double bb_pse, hcal_x, hcal_y;

  T->SetBranchStatus("bb.tr.n",1);
  T->SetBranchStatus("sbs.tr.n",1);
  T->SetBranchStatus("sbs.tr.ph",1);
  T->SetBranchStatus("sbs.tr.th",1);
  T->SetBranchStatus("sbs.tr.y",1);
  T->SetBranchStatus("sbs.tr.x",1);
  T->SetBranchStatus("bb.tr.tg_y",1);
  T->SetBranchStatus("bb.tr.ph",1);
  T->SetBranchStatus("bb.tr.vz",1);
  T->SetBranchStatus("bb.tr.px",1);
  T->SetBranchStatus("bb.tr.py",1);
  T->SetBranchStatus("bb.tr.pz",1);
  T->SetBranchStatus("bb.tr.p",1);
  T->SetBranchStatus("bb.ps.e",1);
  T->SetBranchStatus("sbs.hcal.x",1);
  T->SetBranchStatus("sbs.hcal.y",1);


  T->SetBranchAddress("bb.tr.n",&bb_ntrack);
  T->SetBranchAddress("sbs.tr.n",&sbs_ntrack);
  T->SetBranchAddress("sbs.tr.ph",sbs_ph);
  T->SetBranchAddress("sbs.tr.th",sbs_th);
  T->SetBranchAddress("sbs.tr.y",sbs_y);
  T->SetBranchAddress("sbs.tr.x",sbs_x);
  T->SetBranchAddress("bb.tr.tg_y",bb_tgy);
  T->SetBranchAddress("bb.tr.th",bb_th);
  T->SetBranchAddress("bb.tr.vz",bb_vz);
  T->SetBranchAddress("bb.tr.px",bb_px);
  T->SetBranchAddress("bb.tr.py",bb_py);
  T->SetBranchAddress("bb.tr.pz",bb_pz);
  T->SetBranchAddress("bb.tr.p",bb_p);
  T->SetBranchAddress("bb.ps.e",&bb_pse);
  T->SetBranchAddress("sbs.hcal.x",&hcal_x);
  T->SetBranchAddress("sbs.hcal.y",&hcal_y);

  TH2F *hhBB_SBS_y = new TH2F("hhBB_SBS_y","BB vs SBS y track;SBS track y (m);BB track y (m)",100,-0.3,0.3,100,-0.3,0.3);
  TH2F *hhHCAL_SBS = new TH2F("hhHCAL_SBS","HCal, SBS GEM X corr;SBS GEM #theta (rad);HCal x (m)",100,-0.2,0.1,100,-3,2);
  TH2F *hhBB_SBS = new TH2F("hhBB_SBS","BB, SBS GEM #theta corr;SBS GEM #theta (rad);BB GEM #theta (rad)",100,-0.2,0.1,100,-0.15,0.1);

  int ievent = 0;

  while(T->GetEntry(ievent++)){

    if(sbs_ntrack > 0){
      hhHCAL_SBS->Fill(sbs_th[0],hcal_x);
    }

    if(bb_ntrack > 0 && sbs_ntrack > 0){
      hhBB_SBS_y->Fill(bb_tgy[0],z0_start*tan(sbs_ph[0]) - sbs_y[0]);
      hhBB_SBS->Fill(sbs_th[0],bb_th[0]);
    }
  }
  

  TCanvas *c = new TCanvas("c","",1800,600);
  c->Divide(2,1);
  
  c->cd(1);
  hhHCAL_SBS->Draw("colz");
  c->cd(2);  
  //hhBB_SBS->Draw("colz");
  hhBB_SBS_y->Draw("colz");

}
