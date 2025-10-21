
#include "../../include/gen-ana.h"

void SBS_GEM(){

  double z0_start = 5.0;
  double offset_start = 0.045;

  double ebeam = 4.291;
  double SBS_angle = 34.7*3.1415/180;
  vector<TVector3> GEM_axes; kine::SetHCALaxes(SBS_angle, GEM_axes);
  TVector3 GEM_origin = z0_start*GEM_axes[2];

  TString rootdir = "/lustre19/expphy/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/Rootfiles/pass1/GEN2/H2/SBS100/rootfiles";

  TChain *T = new TChain("T");
  T->Add(rootdir + "/*");

  T->SetBranchStatus("*",0);

  int maxhits = 1000;

  double bb_ntrack, sbs_ntrack;
  double bb_vz[maxhits], bb_tgy[maxhits], bb_px[maxhits], bb_py[maxhits], bb_pz[maxhits], bb_p[maxhits];
  double sbs_ph[maxhits], sbs_y[maxhits], sbs_x[maxhits];
  double bb_pse;

  T->SetBranchStatus("bb.tr.n",1);
  T->SetBranchStatus("sbs.tr.n",1);
  T->SetBranchStatus("sbs.tr.ph",1);
  T->SetBranchStatus("sbs.tr.y",1);
  T->SetBranchStatus("sbs.tr.x",1);
  T->SetBranchStatus("bb.tr.tg_y",1);
  T->SetBranchStatus("bb.tr.vz",1);
  T->SetBranchStatus("bb.tr.px",1);
  T->SetBranchStatus("bb.tr.py",1);
  T->SetBranchStatus("bb.tr.pz",1);
  T->SetBranchStatus("bb.tr.p",1);
  T->SetBranchStatus("bb.ps.e",1);


  T->SetBranchAddress("bb.tr.n",&bb_ntrack);
  T->SetBranchAddress("sbs.tr.n",&sbs_ntrack);
  T->SetBranchAddress("sbs.tr.ph",sbs_ph);
  T->SetBranchAddress("sbs.tr.y",sbs_y);
  T->SetBranchAddress("sbs.tr.x",sbs_x);
  T->SetBranchAddress("bb.tr.tg_y",bb_tgy);
  T->SetBranchAddress("bb.tr.vz",bb_vz);
  T->SetBranchAddress("bb.tr.px",bb_px);
  T->SetBranchAddress("bb.tr.py",bb_py);
  T->SetBranchAddress("bb.tr.pz",bb_pz);
  T->SetBranchAddress("bb.tr.p",bb_p);
  T->SetBranchAddress("bb.ps.e",&bb_pse);

  TH2F *hhcorr = new TH2F("hhcorr","SBS vs BB y track;BB track y (m);SBS track y (m)",100,-0.2,0.2,100,-0.5,0.5);
  TH2F *hhqvec = new TH2F("hhqvec","SBS GEM vs q Vector;#Deltay;#Deltax",100,-0.2,0.2,100,-0.6,0.1);
  TH1F *hqvec = new TH1F("hqvec","Test",100,-0.6,0.1);
  TH1F *hbefore = new TH1F("hbefore","",100,-0.8,0.8);
  TH1F *hafter = new TH1F("hafter","",100,-0.8,0.8);

  int ievent = 0;
  
  while(T->GetEntry(ievent++)){

    if(bb_ntrack == 1 && sbs_ntrack == 1 && abs(bb_vz[0]) < 0.27 && bb_pse > 0.15){
      
      hhcorr->Fill(bb_tgy[0],z0_start*tan(sbs_ph[0]) - sbs_y[0]);

      TVector3 vertex(0,0,bb_vz[0]);
      TLorentzVector Pe(0,0,ebeam,ebeam);   // incoming e-
      TLorentzVector Peprime(bb_px[0],bb_py[0],bb_pz[0],bb_p[0]);   // scattered e-
      double ephi = kine::ephi(Peprime);
      double etheta = kine::etheta(Peprime);
      double pcentral = kine::pcentral(ebeam, etheta, "p");
      double phiN_expect = ephi + constant::pi;

      double nu = Pe.E() - pcentral;
      double pN_expect = kine::pN_expect(nu, "p");
      double thetaN_expect = acos((Pe.E() - pcentral*cos(etheta)) / pN_expect);
      TVector3 pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);

      vector<double> xyGEM_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
      kine::GetxyHCALexpect(vertex, pNhat, GEM_origin, GEM_axes, xyGEM_exp);
    
      hhqvec->Fill(sbs_y[0] - xyGEM_exp[1],sbs_x[0] - xyGEM_exp[0]);
      hqvec->Fill(sbs_x[0] - xyGEM_exp[0]);
      
    }

  }

  TCanvas *c = new TCanvas("c","",1400,600);
  c->Divide(2,1);

  c->cd(1);
  hhcorr->Draw("colz");
  c->cd(2);
  hhqvec->Draw("colz");

  
}
