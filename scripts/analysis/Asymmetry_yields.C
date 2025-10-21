//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified March 4, 2024
//
//
//   The purpose of this script is to calculate the asymmetry
//   given simulated p/n data and real He3 data. This script
//   assumes the data files are already calibrated and analyzed
//   using the script Quasielastic_ana.C
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"

map<int,double> fGoodHel;

// Load database files
void getDB(TString DB_file){

  // Entries should follow this form:
  //{var name, variable pointer, vairable description, 1/0 (mandatory/not mandatory variable)}
  DBparse::DBRequest request[] = {
    {"Good Helicity", &fGoodHel, NULL, "Is the helicity readback good (0/1 = no/yes)", 1}
  };
  
  const int nvar = sizeof(request) / sizeof(request[0]);
  
  DB_load(DB_file,request,nvar);
}

// Input types: GEN2, GEN3, GEN4, GEN4b
void Asymmetry_yields(TString cfg){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TString DB_file = "../../DB/Helicity_quality.csv";
  getDB(DB_file);

  // Load the kinematics and data
  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file);

  // Set variables for cuts
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
  
  // These cuts are for accidental background
  double coin_bg_min = kin_info.coin_time_cut[0] + (1 + 3)*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  double coin_bg_max = kin_info.coin_time_cut[0] + (3 + 3)*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];

  // Set the analyzed tree data
  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_info,1,0);
  analyzed_tree *T_sim = Utilities::LoadAnalyzedRootFiles(kin_info,0,0);

  // Set the distributions
  distribution_fits *dists = new distribution_fits();

  // Set background shape type
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
  
  //dx histograms
  TH1F *hdx_data = new TH1F("hdx_data","",nbins,hxmin,hxmax);
  TH1F *hdx_sim_p = new TH1F("hdx_sim_p","",nbins,hxmin,hxmax);
  TH1F *hdx_sim_n = new TH1F("hdx_sim_n","",nbins,hxmin,hxmax);
  TH1F *hdx_bg_data = new TH1F("hdx_bg_data","",nbins,hxmin,hxmax);
  
  // Cut on the simulated data for QE events from p/n
  // It is also weighted by the "weight" factor from simulation
  TCut CutSimP = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 1) * weight",W2min,W2max,dymin,dymax);
  TCut CutSimN = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 0) * weight",W2min,W2max,dymin,dymax);
  
  // Fill the histograms by drawing
  T_sim->fChain->Draw("dx>>hdx_sim_p",CutSimP);
  T_sim->fChain->Draw("dx>>hdx_sim_n",CutSimN);

  // Asymmetry data for each run
  map<int,analyzed_info*> Asym_runnums;
  
  // Initialize the asymmetry objects for each run
  for(int irun = 0; irun < kin_info.runnums.size(); irun++){
    if(!fGoodHel[kin_info.runnums[irun]]) continue;
    
    analyzed_info *run_info = new analyzed_info();
    Asym_runnums.insert(std::make_pair(kin_info.runnums[irun],run_info));
  }


  int nevent = 0;
  int maxevent = T_data->fChain->GetEntries(); // Used to loop through the tree
  analyzed_info *Asym_total = new analyzed_info(); // Object for total asymmetry

  while(nevent < maxevent){
    T_data->GetEntry(nevent++);  //Get data for one event

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_hel = fGoodHel[T_data->runnum] && (T_data->helicity == -1 || T_data->helicity == 1);
    bool good_W2 = T_data->W2 > W2min && T_data->W2 < W2max;
    bool good_coin_time = T_data->coin_time > coin_min && T_data->coin_time < coin_max;
    bool dy_bg_cut = T_data->dy < dy_bg_min || T_data->dy > dy_bg_max;
    bool good_dy_elas = T_data->dy > dymin && T_data->dy < dymax;
    bool good_dx_elas = T_data->dx > dxmin && T_data->dx < dxmax;
    //////////////////////////////////////////////////////////////////////

    int helicity = T_data->helicity;
    helicity *= -1*T_data->IHWP*IHWP_Flip; 

    if(!good_hel) continue;  //Remove events with bad helicity

    // Fill inelastic helicity data
    if(dy_bg_cut && good_coin_time && good_dx_elas)
      Asym_total->IterateInelasticCount(helicity);

    if(!good_W2) continue;   //Remove events outside of W2 region

    // Use wide dy cut to get bg for fitting
    if(dy_bg_cut && good_coin_time)
      hdx_bg_data->Fill(T_data->dx);  //Fill Bg histogram for fitting later

    if(!good_dy_elas) continue;  //Remove events outside of neutron dy region

    //Cut around bad coin time for bgkd estimation
    if(T_data->coin_time > coin_bg_min && T_data->coin_time < coin_bg_max){
      if(good_dx_elas && good_dy_elas){ // Cut around neutron in dR
	Asym_runnums[T_data->runnum]->IterateBgCount(helicity);
	Asym_total->IterateBgCount(helicity);
      }
    }
    
    if(!good_coin_time) continue;  //Remove events outside of coincidence time

    hdx_data->Fill(T_data->dx);  // Fill data histogram after all QE cuts

    if(good_dx_elas){ //Cut around neutron in dx

      UpdateExpansionCoefficients(T_data);
 
      // Add raw helicity counts
      Asym_total->IterateRawCount(helicity);
      Asym_runnums[T_data->runnum]->IterateRawCount(helicity);
    }
  }

  // Set the distributions used for fitting
  dists->SetDataShape(hdx_data);
  dists->SetPShape(hdx_sim_p);
  dists->SetNShape(hdx_sim_n);
  dists->SetBgShape(hdx_bg_data);
  
  dists->He3_fit_dists();    // Fit our histograms
  
  //Copy all the result histograms
  TH1F *hdx_data_plot = dists->GetDataHist();
  TH1F *hdx_sim_p_plot = dists->GetPHist();
  TH1F *hdx_sim_n_plot = dists->GetNHist();
  TH1F *hdx_bg_plot = dists->GetBgHist();
  TH1F *hdx_total_fit_plot = dists->GetTotalHist();
  
  // Draw historrams
  TCanvas *c = new TCanvas("c","",800,600);
  hdx_data->Draw();
  hdx_data_plot->Draw("hist");
  hdx_sim_p_plot->Draw("same hist");
  hdx_sim_n_plot->Draw("same hist");
  hdx_bg_plot->Draw("same hist");
  
  // Set protona and inelastic yields from the histogram fits
  Asym_total->SetNProton(dists->GetPYield(dxmin, dxmax));
  Asym_total->SetNInelContamination(dists->GetBgYield(dxmin, dxmax));

  double Asym_num = 0;
  double Asym_den = 0;

  // Loop over all asymmetrys run by run
  for (const auto& pair : Asym_runnums) {
    analyzed_info *Asym_results = pair.second;
    
    if(Asym_results->N_raw_p == 0 && Asym_results->N_raw_m == 0) continue;

    /*
    Asym_results->CalcAsymVals();
    
    Asym_num += Asym_results->A_phys / pow(Asym_results->A_phys_stat_err,2);
    Asym_den += 1.0 / pow(Asym_results->A_phys_stat_err,2);
    */
  }
  
  // Calculate the asymmetry from the totals
  Asym_total->CalcAsymVals();

  cout<<Asym_total->f_bg<<" "<<Asym_total->f_p<<" "<<Asym_total->f_in<<endl;
  cout<<Asym_total->A_bg<<" "<<Asym_total->A_p<<" "<<Asym_total->A_in<<endl;


 
}


