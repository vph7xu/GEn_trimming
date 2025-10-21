



void hcal_resolution(TString cfg = "GEN2"){

  gStyle->SetOptStat(0);

  //Read the He3 run and H2 run files
  TFile *H2_file = new TFile("outfiles/QE_test_" + cfg + "_sbs100p_nucleon_p_model2_data.root","read");
  
  TTree *T = (TTree*)H2_file->Get("Tout");


  TH2D *hcorr = new TH2D("hcorr","",100,-2.5,2,100,-0.8,0.8);

  T->Draw("trX:xHCAL>>hcorr","WCut && coinCut","colz");
  
  TF1 *fit = new TF1("fit","pol1",-2.5,0.5);

  hcorr->Fit("fit","R");


  TH1D *hcal_res = new TH1D("hcal_res","",100,-0.3,0.3);
  
  T->Draw(Form("trX - %g*xHCAL - %g>>hcal_res",fit->GetParameter(1),fit->GetParameter(0)),"WCut && coinCut && abs(trX) < 0.1 && abs(trY + 0.1) < 0.002","colz");

  hcal_res->Draw();

  TF1 *fit_gaus = new TF1("fit_gaus","gaus",-0.05,0.05);

  hcal_res->Fit("fit_gaus","R");
  hcal_res->SetTitle(Form("trX - %g*xHCAL - %g",fit->GetParameter(1),fit->GetParameter(0)));
}
