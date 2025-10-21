



void detector_plots(){

  gStyle->SetOptStat(0);

  TString DIR = "/lustre19/expphy/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/Rootfiles/pass1/GEN2/He3/rootfiles/";

  TChain *T = new TChain("T");
  T->Add(DIR + "*2033*");

  TH1F *hps_e = new TH1F("hps_e","Preshower Energy;E (GeV);Counts",100,0,1);
  TH1F *hvz = new TH1F("hvz","Target Vertex Position;Vertex z (m);Counts",100,-0.6,0.6);

  //TCut goodevent = "bb.tr.n==1&&abs(bb.tr.vz[0])<0.27&&sbs.hcal.e>0.025&&bb.ps.e>0.15";
  TCut goodevent = "bb.tr.n==1&&sbs.hcal.e>0.025&&bb.ps.e>0.2";

  T->Draw("bb.ps.e>>hps_e");
  T->Draw("bb.tr.vz[0]>>hvz");


  TCanvas *c = new TCanvas("c","",800,600);

  hps_e->Draw();

  TLine *line = new TLine(0.2,0,0.2,hps_e->GetMaximum());
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("same");

  TCanvas *c2 = new TCanvas("c2","",800,600);

  hvz->Draw();

  TLine *line2 = new TLine(-0.27,0,-0.27,8600);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(2);
  line2->Draw("same");

  TLine *line3 = new TLine(0.27,0,0.27,8600);
  line3->SetLineColor(kRed);
  line3->SetLineWidth(2);
  line3->Draw("same");

 
  TString kin_name[3] = {"GEN2","GEN3","GEN4"};
  TString title[3] = {"Q^{2} = 2.9 GeV","Q^{2} = 6.6 GeV","Q^{2} = 9.7 GeV"};
  TH2F *h2[3];
  TH1F *h1[3];

  

  for(int ikin=0; ikin < 3; ikin++){
    
    TFile *file = new TFile("../outfiles/QE_data_" + kin_name[ikin] + "_sbs100p_nucleon_np_model2.root","read");
    TTree *T = (TTree*)file->Get("Tout");

    h2[ikin] = new TH2F("h2_" + kin_name[ikin],title[ikin] + ";W^{2} (GeV);#Deltax (m)",100,-1,5,100,-4,1);
    h1[ikin] = new TH1F("h1_" + kin_name[ikin],title[ikin] + ";W^{2} (GeV)",100,-1,5);

    T->Draw("dx:W2>>h2_" + kin_name[ikin]);
    T->Draw("W2>>h1_" + kin_name[ikin]);
    h2[ikin]->GetYaxis()->SetTitleOffset(1.2);
    
  }


  TCanvas *c5 = new TCanvas("c5","",1500,600);
  c5->Divide(3,1);

  c5->cd(1);
  h2[0]->Draw("colz");
  c5->cd(2);
  h2[1]->Draw("colz");
  c5->cd(3);
  h2[2]->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","",800,600);
  h1[0]->Draw();

}
