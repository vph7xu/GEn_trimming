
const int nQ2 = 3;
const int nlayers = 5;

double Q2_list[nQ2] = {3.0,6.8,9.8};
TString kin_list[nQ2] = {"GEN2","GEN3","GEN4"};
int runs[nQ2] = {2164,2566,3854};

int NStripsU[nlayers] = {3840, 3840, 3840, 3840, 5120};
int NStripsV[nlayers] = {3840, 3840, 3840, 3840, 1536*4};


double get_eff(TFile *file, int layer){
  TH2F *hdid = (TH2F*)file->Get(Form("hdidhit_xy_bb_gem_layer%i",layer));
  TH2F *hshould = (TH2F*)file->Get(Form("hshouldhit_xy_bb_gem_layer%i",layer));

  return hdid->GetEntries()*1.0/hshould->GetEntries();
}

void occupancy_plots(){

  double U_occu[nQ2][nlayers];
  double U_err[nQ2][nlayers];
  double efficiency[nQ2][nlayers];


  TString rootdir = "/cache/halla/sbs/prod/GEnII/pass1/";

  for(int iQ2=0; iQ2 < nQ2; iQ2++){
    
    TChain *T = new TChain("T");
    T->Add(Form(rootdir + kin_list[iQ2] + "/He3/rootfiles/*%i*",runs[iQ2]));
    
    TFile *file = new TFile(Form(rootdir + kin_list[iQ2] + "/He3/rootfiles/e1209016_fullreplay_%i_stream0_2_seg0_0.root",runs[iQ2]),"read");

    for(int ilayer = 0; ilayer < nlayers; ilayer++){
      TH1F *hist = new TH1F("hist","",100,0,1000);
    
      T->Draw(Form("bb.gem.nstripsu_layer[%i]>>hist",ilayer));
    
      TF1 *fit = new TF1("fit","gaus");
	  
      hist->Fit(fit,"q0");
      
      U_occu[iQ2][ilayer] = fit->GetParameter(1) / NStripsU[ilayer];
      U_err[iQ2][ilayer] = fit->GetParameter(2) / NStripsU[ilayer];
      
      hist->Delete();

      efficiency[iQ2][ilayer] = get_eff(file,ilayer);
    }
  }
  
  TGraphErrors *g_occu[nQ2];
  TGraph *g_eff[nQ2];
  
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TCanvas *c2 = new TCanvas("c2","",800,600);
  TLegend *legend = new TLegend(0.60,0.70,0.9,0.9);
  int icolor = 0;

  double layer_list[nlayers] = {1,2,3,4,5};

  for(int iQ2 = 0; iQ2 < nQ2; iQ2++){
    
    icolor++;

    g_occu[iQ2] = new TGraphErrors(nlayers,layer_list,U_occu[iQ2],0,U_err[iQ2]);
    g_occu[iQ2]->SetTitle("GEN GEM Layer Occupancy;Layer # ;Occupancy");
    g_occu[iQ2]->SetMarkerStyle(8);
    g_occu[iQ2]->SetMarkerColor(icolor);
    g_occu[iQ2]->SetLineColor(icolor);
    legend->AddEntry(g_occu[iQ2],Form("Q^{2} = %g",Q2_list[iQ2]),"p");

    g_eff[iQ2] = new TGraph(nlayers,layer_list,efficiency[iQ2]);
    g_eff[iQ2]->SetTitle("GEN GEM Layer Efficiency;Layer # ;Efficiency");
    g_eff[iQ2]->SetMarkerStyle(8);
    g_eff[iQ2]->SetMarkerColor(icolor);
    g_eff[iQ2]->SetLineColor(icolor);

    if(iQ2 == 0){
      c1->cd();
      g_occu[iQ2]->Draw("AP");
      legend->Draw("same");

      c2->cd();
      g_eff[iQ2]->Draw("AP");
      legend->Draw("same");
    }
    else{
      c1->cd();
      g_occu[iQ2]->Draw("P");

      c2->cd();
      g_eff[iQ2]->Draw("P");
    }    

    g_occu[iQ2]->GetYaxis()->SetRangeUser(0,0.2);
    g_eff[iQ2]->GetYaxis()->SetRangeUser(0.4,1.1);
    
  }


}
