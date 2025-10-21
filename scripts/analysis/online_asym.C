//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is to check the asymmetry 
//   using pion events. It will plot the asymmetry as a number
//   of events.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#include "../../include/gen-ana.h"


void online_asym(const char *configfilename){

  gStyle->SetOptStat(0);

  int maxevent = 3000000;
  
  const int npoints = 10;

  double maxevent_test[npoints];

  for(int i=0;i < npoints; i++) maxevent_test[i] = 1 + i*1;

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
  TChain *C = new TChain("T");
  if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();
  for (int i=0; i<nruns; i++) {
    std::string rfname = rootfile_dir + Form("/*%d*",runnums[i]);
    C->Add(rfname.c_str());
  }

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("GEN_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);

  int IHWP_Flip = jmgr->GetValueFromKey<int>("IHWP_Flip");

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  int model = jmgr->GetValueFromKey<int>("model");
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromKey_str("Ntype");


  C->SetBranchStatus("*",0);

  //Beam helicity variables
  double helicity;
  setrootvar::setbranch(C,"scalhel","hel",&helicity);

  //IHWP State
  double IHWP;
  setrootvar::setbranch(C,"IGL1I00OD16_16","",&IHWP); 

  // track var
  int maxNtr=1000;
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr],xTr[maxNtr],yTr[maxNtr],thTr[maxNtr],phTr[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&xTr,&yTr,&thTr,&phTr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  // hcal clus var
  double eHCAL, xHCAL, yHCAL, rblkHCAL, cblkHCAL, idblkHCAL,tdctimeHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","tdctimeblk"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&tdctimeHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  // bbcal clus var
  double eSH,atimeSH;
  std::vector<std::string> bbcalclvar = {"e","atimeblk"}; 
  std::vector<void*> bbcalclvar_mem = {&eSH,&atimeSH}; 
  setrootvar::setbranch(C,"bb.sh",bbcalclvar,bbcalclvar_mem);
  double ePS; setrootvar::setbranch(C,"bb.ps","e",&ePS);

  // hodoscope
  const int maxClus = 1000;
  double hodo_time[maxClus]; setrootvar::setbranch(C,"bb.hodotdc.clus.bar.tdc","meantime",&hodo_time);
  int nhodo_clus; setrootvar::setbranch(C,"Ndata.bb.hodotdc.clus.bar.tdc","meantime",&nhodo_clus);


  TH1D *h1 = new TH1D("h1","Photon Energy Spectrum from Deteced Pions;E_{#gamma};Counts",150,0,15);
  TH1D *hps = new TH1D("hps","BB Preshower Energy;E_{ps} (GeV);Counts",150,0,1);
  TH2D *hcorr = new TH2D("hcorr","BB #pi- correlation to HCAL;BB tr #phi (rad);HCAL x (m)",150,-0.6,0.6,150,-3,2);
  //TH2D *hcorr = new TH2D("hcorr","BB #pi- correlation to HCAL;BB tr #phi (rad);HCAL x (m)",150,1,2,150,-3,2);

  

  int Yp_e = 0;
  int Ym_e = 0;
  int Yp_pim = 0;
  int Ym_pim = 0;
  int Yp_pi0 = 0;
  int Ym_pi0 = 0;

  int Yp_pi_tot[npoints] = {0};
  int Ym_pi_tot[npoints] = {0};

  for(int i=0; i < npoints; i++){
    
    maxevent = maxevent_test[i]*1000000;
  
    int nevent = 0;  

    Yp_e = 0;
    Ym_e = 0;
    Yp_pim = 0;
    Ym_pim = 0;
    Yp_pi0 = 0;
    Ym_pi0 = 0;
  while (C->GetEntry(nevent++) && nevent < maxevent) {
    
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent*100.0/maxevent << "% \r";
    std::cout.flush();

    double coin_time = hodo_time[0] - tdctimeHCAL;
    
    //bool goodevent = ntrack==1 && abs(vz[0]) < 0.27 && eHCAL > 0.025 && coin_time > 3.5 && coin_time < 16.3;
    bool goodevent = ntrack==1 && abs(vz[0]) < 0.27 && eHCAL > 0.025;

    //if(!goodevent) continue;

    double ebeam = sbsconf.GetEbeam();
    TLorentzVector Ppiprime(px[0],   // scattered pi-
			   py[0],
			   pz[0],
			   p[0]);

    TLorentzVector Pn(0,0,0,constant::Mn);    //Rest neutron

    double E_gamma = (2*constant::Mn*Ppiprime.E() + pow(constant::Mp,2) - pow(constant::Mpim,2) - pow(constant::Mn,2)) / (2*(constant::Mn + Ppiprime.Px() - Ppiprime.E()));

    TLorentzVector Pgamma(0,0,E_gamma,E_gamma);    //incoming photon

    TLorentzVector Ppprime = Ppiprime - Pgamma - Pn;

    h1->Fill(E_gamma);
    hps->Fill(ePS);


    if(IHWP == 1) helicity *= -1*IHWP_Flip;
    else if(IHWP == -1) helicity *= 1*IHWP_Flip;
    else continue;


    bool piCut = ePS < 0.2 && ntrack==1 && abs(vz[0]) < 0.27 && eHCAL > 0.025;
    bool eCut = ePS > 0.2 && ntrack==1 && abs(vz[0]) < 0.27 && eHCAL > 0.025;
    bool pi0Cut = ePS < 0.2 && ePS > 0.04 && ntrack==0 && eHCAL > 0.025;
    
    if(eCut){
      if(helicity == 1) Yp_e++;
      if(helicity == -1) Ym_e++;
    }

    if(piCut){
      if(E_gamma > 0.75*ebeam && E_gamma < 1.25*ebeam) hcorr->Fill(atan2(py[0],px[0]),xHCAL);
      if(helicity == 1){
	Yp_pim++;
	Yp_pi_tot[i]++;
      }
      if(helicity == -1){
	Ym_pim++;
	Ym_pi_tot[i]++;
      }
    }
      
    if(pi0Cut){
      if(helicity == 1){
	Yp_pi0++;
	Yp_pi_tot[i]++;
      }
      if(helicity == -1){
	Ym_pi0++;
	Ym_pi_tot[i]++;
      }
    }
  }
  }
 

  //h1->Draw();
  //hcorr->Draw("colz");
  //hps->Draw();

  double A_e = (Yp_e - Ym_e)*1.0/(Yp_e + Ym_e);
  double A_e_err = sqrt((1 - A_e*A_e)/(Yp_e + Ym_e));

  double A_pim = (Yp_pim - Ym_pim)*1.0/(Yp_pim + Ym_pim);
  double A_pim_err = sqrt((1 - A_pim*A_pim)/(Yp_pim + Ym_pim));

  double A_pi0 = (Yp_pi0 - Ym_pi0)*1.0/(Yp_pi0 + Ym_pi0);
  double A_pi0_err = sqrt((1 - A_pi0*A_pi0)/(Yp_pi0 + Ym_pi0));

  double Yp_pi = Yp_pim + Yp_pi0;
  double Ym_pi = Ym_pim + Ym_pi0;

  double A_pi = (Yp_pi - Ym_pi)*1.0/(Yp_pi + Ym_pi);
  double A_pi_err = sqrt((1 - A_pi*A_pi)/(Yp_pi + Ym_pi));

  double A_pi_tot[npoints], A_pi_tot_err[npoints];
  
  for(int i=0; i<npoints; i++){
    A_pi_tot[i] = (Yp_pi_tot[i] - Ym_pi_tot[i])*1.0/(Yp_pi_tot[i] + Ym_pi_tot[i]);
    A_pi_tot_err[i] = sqrt((1 - A_pi_tot[i]*A_pi_tot[i])/(Yp_pi_tot[i] + Ym_pi_tot[i]));
    A_pi_tot[i] *= 100;
    A_pi_tot_err[i] *= 100;
}

  TGraphErrors *gA = new TGraphErrors(npoints,maxevent_test,A_pi_tot,0,A_pi_tot_err);
  
  TCanvas *c = new TCanvas("c","",800,600);
  gA->Draw("AP");
  gA->SetTitle("Pion Statistics Needed for Asymmetry;Replayed Events (Millions);Asymmetry (%)");
  gA->SetMarkerStyle(8);

  cout<<"Total electron events = "<<Yp_e + Ym_e<<endl;
  cout<<"Y+ = "<<Yp_e<<endl;
  cout<<"Y- = "<<Ym_e<<endl;
  cout<<"electron A = "<<A_e<<" +/- "<<A_e_err<<endl;
  cout<<"Total pi- events = "<<Yp_pim + Ym_pim<<endl;
  cout<<"Y+ = "<<Yp_pim<<endl;
  cout<<"Y- = "<<Ym_pim<<endl;
  cout<<"pi- A = "<<A_pim<<" +/- "<<A_pim_err<<endl;
  cout<<"Total pi0 events = "<<Yp_pi0 + Ym_pi0<<endl;
  cout<<"Y+ = "<<Yp_pi0<<endl;
  cout<<"Y- = "<<Ym_pi0<<endl;
  cout<<"pi0 A = "<<A_pi0<<" +/- "<<A_pi0_err<<endl;
  cout<<"Total pi events = "<<Yp_pi + Ym_pi<<endl;
  cout<<"Y+ = "<<Yp_pi<<endl;
  cout<<"Y- = "<<Ym_pi<<endl;
  cout<<"Total pi A = "<<A_pi<<" +/- "<<A_pi_err<<endl;




}
