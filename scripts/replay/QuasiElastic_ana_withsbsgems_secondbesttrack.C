//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Provakar Datta
//   Modified by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is to take a configuraiton
//   file for some SBS experiment and to produced some 
//   analyzed output with cuts on good elastic events.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../include/gen-ana.h"


DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){
  
  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
     {"He3 Polarization","He3 target polarization", 1},
     {"Beam Polarization","Beam Polarization values",1},
     {"Helicity Quality","Helicity readback good? (0/1 = bad/good)",1},
     {"Moller Quality","Moller measurements known? (0/1 = no/yes)",1},
     {"Field Measurement","Magnetic field direction",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}


int QuasiElastic_ana_withsbsgems(const std::string configfilename, std::string filebase="../outfiles/pass2_0731/QE_data")
{

  string configdir = "../../config/";

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(configdir + configfilename,1);

  getDB(kin_info.conf);

  // parsing trees
  TChain *C = LoadRawRootFiles(kin_info, 1);

  // seting up the desired SBS configuration
  TString conf = kin_info.conf;
  int sbsmag = kin_info.sbsmag;
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  int model = kin_info.model;
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }

  // choosing nucleon type 
  TString Ntype = kin_info.Ntype;

  double GEMpitch = 10.0*TMath::DegToRad();
  TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
  TVector3 GEMyaxis(0,1,0);
  TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();

  //sbs gem axes
  

  // setting up ROOT tree branch addresses ---------------------------------------

  // setting up global cuts
  TCut globalcut = kin_info.globalcut;
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // beam energy - Probably we should take an average over 100 events
  // double HALLA_p;
  // setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

  // Some CODA event information
  double evtime;   setrootvar::setbranch(C,"g","evtime",&evtime);
  
  // bbcal sh clus var
  double eSH,xSH,ySH,atimeSH;
  std::vector<std::string> bbcalclvar = {"e","x","y","atimeblk"}; 
  std::vector<void*> bbcalclvar_mem = {&eSH,&xSH,&ySH,&atimeSH}; 
  setrootvar::setbranch(C,"bb.sh",bbcalclvar,bbcalclvar_mem);
  
  // bbcal ps clus var
  double ePS, xPS;
  std::vector<std::string> bbcalpsclvar = {"e","x"}; 
  std::vector<void*> bbcalpsclvar_mem = {&ePS,&xPS}; 
  setrootvar::setbranch(C,"bb.ps",bbcalpsclvar,bbcalpsclvar_mem);
  
  int maxhcal = 100;

  // hcal clus var
  double eHCAL[maxhcal],xHCAL[maxhcal], yHCAL[maxhcal], rblkHCAL[maxhcal], cblkHCAL[maxhcal], idblkHCAL[maxhcal],tdctimeHCAL[maxhcal],atimeHCAL[maxhcal],nclusHCAL,nblkHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","tdctimeblk","atimeblk","nclus","nblk"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&tdctimeHCAL,&atimeHCAL,&nclusHCAL,&nblkHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  //hcal all clus var
  double clusEHCAL[maxhcal], clustdctimeHCAL[maxhcal], blkid[maxhcal], nblk[maxhcal];
  std::vector<std::string> hcalclmemvar = {"e","tdctime","id","nblk"};
  std::vector<void*> hcalclmemvar_mem = {&clusEHCAL,&clustdctimeHCAL,&blkid,&nblk};
  setrootvar::setbranch(C, "sbs.hcal.clus", hcalclmemvar, hcalclmemvar_mem);

  //hcal primary clus var
  double clusmemEHCAL[maxhcal], clusmematimeHCAL[maxhcal], clusmemtdctimeHCAL[maxhcal], clusmemblkid[maxhcal];
  std::vector<std::string> hcalclusmemvar = {"e","atime","tdctime","id"};
  std::vector<void*> hcalclusmemvar_mem = {&clusmemEHCAL, &clusmematimeHCAL, &clusmemtdctimeHCAL, &clusmemblkid};
  setrootvar::setbranch(C, "sbs.hcal.clus_blk", hcalclusmemvar, hcalclusmemvar_mem);

  // grinch var
  double grinch_track,grinch_clus_size;
  std::vector<std::string> grinchvar = {"trackindex","size"};
  std::vector<void*> grinchvar_mem = {&grinch_track,&grinch_clus_size};
  setrootvar::setbranch(C, "bb.grinch_tdc.clus", grinchvar, grinchvar_mem);
  

  // hodoscope
  const int maxClus = 1000;
  double hodo_time[maxClus]; 
  int nhodo_clus; 
  
  setrootvar::setbranch(C,"bb.hodotdc.clus.bar.tdc","meantime",&hodo_time);
  setrootvar::setbranch(C,"Ndata.bb.hodotdc.clus.bar.tdc","meantime",&nhodo_clus);

  //sbsgem
  const int maxTracks = 1000;
  double sbs_gem_nhits[maxTracks];
  int nsbs_gem_nhits;

  double bb_gem_nhits[maxTracks];
  int nbb_gem_nhits;

  setrootvar::setbranch(C,"sbs.gem.track","nhits",&sbs_gem_nhits);
  setrootvar::setbranch(C,"bb.gem.track","nhits",&bb_gem_nhits);  
  //setrootvar::setbranch(C,"Ndata.sbs.gem.track","nhits",&nsbs_gem_nhits);


  // track var
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr],xTr[maxNtr],yTr[maxNtr],thTr[maxNtr],phTr[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr];
  double xfp[maxNtr],yfp[maxNtr],thfp[maxNtr],phfp[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph","r_x","r_y","r_th","r_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,&xfp,&yfp,&thfp,&phfp};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  //sbs track var
  double ntrack_sbs, p_sbs[maxNtr],px_sbs[maxNtr],py_sbs[maxNtr],pz_sbs[maxNtr],xTr_sbs[maxNtr],yTr_sbs[maxNtr],thTr_sbs[maxNtr],phTr_sbs[maxNtr];
  double vx_sbs[maxNtr],vy_sbs[maxNtr],vz_sbs[maxNtr];
  double ytgt_sbs[maxNtr],thtgt_sbs[maxNtr],phtgt_sbs[maxNtr];
  double xfp_sbs[maxNtr],yfp_sbs[maxNtr],thfp_sbs[maxNtr],phfp_sbs[maxNtr];
  std::vector<std::string> sbstrvar ={"n","p","px","py","pz","vx","vy","vz","tg_y","tg_th","tg_ph","r_x","r_y","r_th","r_ph","x","y","th","ph"};
  std::vector<void*> sbstrvar_mem = {&ntrack_sbs,&p_sbs,&px_sbs,&py_sbs,&pz_sbs,&vx_sbs,&vy_sbs,&vz_sbs,&ytgt_sbs,&thtgt_sbs,&phtgt_sbs,&xfp_sbs,&yfp_sbs,&thfp_sbs,&phfp_sbs,&xTr_sbs,&yTr_sbs,&thTr_sbs,&phTr_sbs};
  setrootvar::setbranch(C,"sbs.tr",sbstrvar,sbstrvar_mem);

  // tdctrig variable (N/A for simulation)
  int tdcElemN;
  double tdcTrig[maxNtr], tdcElem[maxNtr];
  
  std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
  std::vector<void*> tdcvar_mem = {&tdcElem,&tdcElemN,&tdcTrig};
  setrootvar::setbranch(C,"bb.tdctrig",tdcvar,tdcvar_mem,1);
  

  //Beam helicity variables
  double helicity;
  setrootvar::setbranch(C,"scalhel","hel",&helicity);

  //IHWP State
  double IHWP;
  setrootvar::setbranch(C,"IGL1I00OD16_16","",&IHWP); 

  //trigbits
  double trigbits;
  setrootvar::setbranch(C,"g","trigbits",&trigbits); 

  //BPMA
  double BPMAx, BPMAy;
  std::vector<std::string> BPMAvar = {"x","y"};
  std::vector<void*> BPMAvar_mem = {&BPMAx,&BPMAy};
  setrootvar::setbranch(C,"Lrb.BPMA",BPMAvar,BPMAvar_mem);

  //Raster
  double rawcurx, rawcury;
  std::vector<std::string> Rastervar = {"x","y"};
  std::vector<void*> Rastervar_mem = {&rawcurx,&rawcury};
  setrootvar::setbranch(C,"Lrb.Raster.rawcur",Rastervar,Rastervar_mem);

  //Raster 2
  double rawcur2x, rawcur2y;
  std::vector<std::string> Raster2var = {"x","y"};
  std::vector<void*> Raster2var_mem = {&rawcur2x,&rawcur2y};
  setrootvar::setbranch(C,"Lrb.Raster2.rawcur",Raster2var,Raster2var_mem);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);
  C->SetBranchStatus("sbs.hcal.nclus", 1);
  
  // defining the outputfile
  TString outFile = Form("%s_" + sbsconf.GetSBSconf() + "_sbs%dp_nucleon_%s_model%d_sbstrackingon.root", 
			 filebase.c_str(),  sbsconf.GetSBSmag(), Ntype.Data(), model);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining histograms
  TH1F *h_W = Utilities::TH1FhW("h_W");
  TH1F *h_W_cut = Utilities::TH1FhW("h_W_cut");
  TH1F *h_W_acut = Utilities::TH1FhW("h_W_acut");
  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  
  TH1F *h_Q2 = Utilities::TH1FhQ2("h_Q2", conf);
  vector<double> hdx_lim = kin_info.hdx_lim;
  vector<double> hdy_lim = kin_info.hdy_lim;
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_coin_time = new TH1F("h_coin_time", "Coincidence time (ns)", 200, 380, 660);

  TH2F *h2_rcHCAL = Utilities::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = Utilities::TH2FdxdyHCAL("h2_dxdyHCAL");
  TH2F *h2_xyHCAL_p = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_p");
  TH2F *h2_xyHCAL_n = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_n");

  // Defining interesting ROOT tree branches 
  TTree *Tout = new TTree("Tout", "");
  int T_runnum;         Tout->Branch("runnum", &T_runnum, "runnum/I");
  TDatime T_datetime;   Tout->Branch("datetime", "TDatime", &T_datetime);
  //cuts
  bool WCut;            Tout->Branch("WCut", &WCut, "WCut/B");
  bool pCut;            Tout->Branch("pCut", &pCut, "pCut/B");
  bool nCut;            Tout->Branch("nCut", &nCut, "nCut/B");
  bool fiduCut;         Tout->Branch("fiduCut", &fiduCut, "fiduCut/B");
  bool coinCut;         Tout->Branch("coinCut", &coinCut, "coinCut/B");
  //
  double T_ebeam;       Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
  //kine
  double T_nu;          Tout->Branch("nu", &T_nu, "nu/D");
  double T_Q2;          Tout->Branch("Q2", &T_Q2, "Q2/D");
  double T_W2;          Tout->Branch("W2", &T_W2, "W2/D");
  double T_dpel;        Tout->Branch("dpel", &T_dpel, "dpel/D");
  double T_ephi;        Tout->Branch("ephi", &T_ephi, "ephi/D");
  double T_etheta;      Tout->Branch("etheta", &T_etheta, "etheta/D");
  double T_pcentral;    Tout->Branch("pcentral", &T_pcentral, "pcentral/D");
  double T_pN_expect; 	Tout->Branch("pN_expect", &T_pN_expect, "pN_expect/D");
  
  double T_ptheta_cal;	Tout->Branch("ptheta_cal",&T_ptheta_cal,"ptheta_cal/D");
  double T_ptheta; 	Tout->Branch("ptheta",&T_ptheta,"ptheta/D");
  double T_pphi_cal;  Tout->Branch("pphi_cal",&T_pphi_cal,"pphi_cal/D");
  double T_pphi;      Tout->Branch("pphi",&T_pphi,"pphi/D");
  
  //track
  double T_ntrack;	Tout->Branch("ntrack", &T_ntrack, "ntrack/D");
  double T_ntrack_hits;	Tout->Branch("ntrack_hits", &T_ntrack_hits, "ntrack_hits/D");
  double T_vz;          Tout->Branch("vz", &T_vz, "vz/D");
  double T_vx;          Tout->Branch("vx", &T_vx, "vx/D");
  double T_vy;          Tout->Branch("vy", &T_vy, "vy/D");
  double T_xtgt;        Tout->Branch("xtgt", &T_xtgt, "xtgt/D");
  double T_ytgt;        Tout->Branch("ytgt", &T_ytgt, "ytgt/D");
  double T_thtgt;       Tout->Branch("thtgt", &T_thtgt, "thtgt/D");
  double T_phtgt;       Tout->Branch("phtgt", &T_phtgt, "phtgt/D");
  double T_thetabend;   Tout->Branch("thetabend", &T_thetabend, "thetabend/D");
  double T_xfp;         Tout->Branch("xfp", &T_xfp, "xfp/D");
  double T_yfp;         Tout->Branch("yfp", &T_yfp, "yfp/D");
  double T_thfp;        Tout->Branch("thfp", &T_thfp, "thfp/D");
  double T_phfp;        Tout->Branch("phfp", &T_phfp, "phfp/D");  
  double T_trP;         Tout->Branch("trP", &T_trP, "trP/D");
  double T_trPx;         Tout->Branch("trPx", &T_trPx, "trPx/D");
  double T_trPy;         Tout->Branch("trPy", &T_trPy, "trPy/D");
  double T_trPz;         Tout->Branch("trPz", &T_trPz, "trPz/D");

  //sbs track 
  double T_ntrack_sbs;       Tout->Branch("ntrack_sbs", &T_ntrack_sbs, "ntrack_sbs/D");
  double T_ntrack_hits_sbs;	Tout->Branch("ntrack_hits_sbs", &T_ntrack_hits_sbs, "ntrack_hits_sbs/D");
  double T_vz_sbs;          Tout->Branch("vz_sbs", &T_vz_sbs, "vz_sbs/D");
  double T_vx_sbs;          Tout->Branch("vx_sbs", &T_vx_sbs, "vx_sbs/D");
  double T_vy_sbs;          Tout->Branch("vy_sbs", &T_vy_sbs, "vy_sbs/D");
  //double T_xtgt_sbs;        Tout->Branch("xtgt_sbs", &T_xtgt_sbs, "xtgt_sbs/D");
  double T_ytgt_sbs;        Tout->Branch("ytgt_sbs", &T_ytgt_sbs, "ytgt_sbs/D");
  double T_thtgt_sbs;       Tout->Branch("thtgt_sbs", &T_thtgt_sbs, "thtgt_sbs/D");
  double T_phtgt_sbs;       Tout->Branch("phtgt_sbs", &T_phtgt_sbs, "phtgt_sbs/D");
  //double T_thetabend_sbs;   Tout->Branch("thetabend_sbs", &T_thetabend_sbs, "thetabend_sbs/D");
  double T_xfp_sbs;         Tout->Branch("xfp_sbs", &T_xfp_sbs, "xfp_sbs/D");
  double T_yfp_sbs;         Tout->Branch("yfp_sbs", &T_yfp_sbs, "yfp_sbs/D");
  double T_thfp_sbs;        Tout->Branch("thfp_sbs", &T_thfp_sbs, "thfp_sbs/D");
  double T_phfp_sbs;        Tout->Branch("phfp_sbs", &T_phfp_sbs, "phfp_sbs/D");
  double T_trP_sbs;         Tout->Branch("trP_sbs", &T_trP_sbs, "trP_sbs/D");
  double T_trPx_sbs;         Tout->Branch("trPx_sbs", &T_trPx_sbs, "trPx_sbs/D");
  double T_trPy_sbs;         Tout->Branch("trPy_sbs", &T_trPy_sbs, "trPy_sbs/D");
  double T_trPz_sbs;         Tout->Branch("trPz_sbs", &T_trPz_sbs, "trPz_sbs/D");
  double T_trx_sbs;          Tout->Branch("trx_sbs",&T_trx_sbs,"trx_sbs/D");
  double T_try_sbs;	     Tout->Branch("try_sbs",&T_try_sbs,"try_sbs/D");
  double T_trth_sbs;	     Tout->Branch("trth_sbs",&T_trth_sbs,"trth_sbs/D");
  double T_trph_sbs;	     Tout->Branch("trph_sbs",&T_trph_sbs,"trph_sbs/D");

  //BBCAL
  double T_ePS;         Tout->Branch("ePS", &T_ePS, "ePS/D"); 
  double T_xPS;         Tout->Branch("xPS", &T_xPS, "xPS/D"); 
  double T_eSH;         Tout->Branch("eSH", &T_eSH, "eSH/D"); 
  double T_xSH;         Tout->Branch("xSH", &T_xSH, "xSH/D"); 
  double T_ySH;         Tout->Branch("ySH", &T_ySH, "ySH/D"); 

  //HCAL
  double T_eHCAL;       Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); 
  //double T_e_cHCAL;	Tout->Branch("e_cHCAL",&T_e_cHCAL,"e_cHCAL/D");
  double T_xHCAL;       Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); 
  double T_yHCAL;       Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D"); 
  double T_xHCAL_exp;   Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D"); 
  double T_yHCAL_exp;   Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D"); 
  double T_dx;          Tout->Branch("dx", &T_dx, "dx/D"); 
  double T_dy;          Tout->Branch("dy", &T_dy, "dy/D");
  double T_theta_pq;    Tout->Branch("theta_pq", &T_theta_pq, "theta_pq/D");
  double T_nclus_HCAL;	Tout->Branch("nclus_HCAL", &T_nclus_HCAL, "nclus_HCAL/D");
  double T_nblk_HCAL;	Tout->Branch("nblk_HCAL", &T_nblk_HCAL, "nblk_HCAL/D");

  //HCAL all clus
  //double T_hcal_clus_E[1000]; Tout->Branch("hcal_clus_e", &T_hcal_clus_E, "hcal_clus_e[1000]/D");
  //double T_hcal_clus_atime[1000]; Tout->Branch("hcal_clus_atime", &T_hcal_clus_atime, "hcal_clus_atime[1000]/D");
  //double T_hcal_clus_tdctime[1000]; Tout->Branch("hcal_clus_tdctime", &T_hcal_clus_tdctime, "hcal_clus_tdctime[1000]/D");
  //double T_hcal_clus_id[1000]; Tout->Branch("hcal_clus_id", &T_hcal_clus_id, "hcal_clus_id[1000]/D");
  //double T_hcal_clus_nblk[1000]; Tout->Branch("hcal_clus_nblk", &T_hcal_clus_nblk, "hcal_clus_nblk[1000]/D");

  //HCAL clus mem
  //double T_hcal_clus_mem_E[1000]; Tout->Branch("hcal_clus_mem_e", &T_hcal_clus_mem_E, "hcal_clus_mem_e[1000]/D");
  //double T_hcal_clus_mem_atime[1000]; Tout->Branch("hcal_clus_mem_atime", &T_hcal_clus_mem_atime, "hcal_clus_mem_atime[1000]/D");
  //double T_hcal_clus_mem_tdctime[1000]; Tout->Branch("hcal_clus_mem_tdctime", &T_hcal_clus_mem_tdctime, "hcal_clus_mem_tdctime[1000]/D");
  //double T_hcal_clus_mem_id[1000]; Tout->Branch("hcal_clus_mem_id", &T_hcal_clus_mem_id, "hcal_clus_mem_id[1000]/D");
  //double T_hcal_clus_mem_nblk[1000]; Tout->Branch("hcal_clus_mem_nblk", &T_hcal_clus_mem_nblk, "hcal_clus_mem_nblk[1000]/D");

  //GRINCH
  double T_grinch_track;    Tout->Branch("grinch_track", &T_grinch_track, "grinch_track/D");
  double T_grinch_clus_size;    Tout->Branch("grinch_clus_size", &T_grinch_clus_size, "grinch_clus_size/D");
  
  //Timing Information
  double T_coin_time;           Tout->Branch("coin_time", &T_coin_time, "coin_time/D");
  double T_hcal_time;            Tout->Branch("hcal_time", &T_hcal_time, "hcal_time/D"); 
  double T_bbcal_time;           Tout->Branch("bbcal_time", &T_bbcal_time, "bbcal_time/D");
  int T_nhodo_clus;              Tout->Branch("nhodo_clus", &T_nhodo_clus, "nhodo_clus/I");  
  double T_hodo_time[maxClus];   Tout->Branch("hodo_time", &T_hodo_time, "hodo_time[nhodo_clus]/D");

  //BPM and Raster information
  double T_BPMAx;     Tout->Branch("BPMAx", &T_BPMAx, "BPMAx/D");
  double T_BPMAy;     Tout->Branch("BPMAy", &T_BPMAy, "BPMAy/D");
  double T_rawcurx;  Tout->Branch("Rasterx", &T_rawcurx, "Rasterx/D");
  double T_rawcury;  Tout->Branch("Rastery", &T_rawcury, "Rastery/D");
  double T_rawcur2x;  Tout->Branch("Raster2x", &T_rawcur2x, "Raster2x/D");
  double T_rawcur2y;  Tout->Branch("Raster2y", &T_rawcur2y, "Raster2y/D");

  //trigbits
  double T_trigbits; Tout->Branch("trigbits", &T_trigbits, "trigbits/D");

  //Beam/Target information
  int T_helicity;        Tout->Branch("helicity", &T_helicity, "helicity/I");
  int T_IHWP;            Tout->Branch("IHWP", &T_IHWP, "IHWP/I");
  double T_He3Pol;       Tout->Branch("He3Pol", &T_He3Pol, "He3Pol/D");
  double T_err_He3Pol;       Tout->Branch("err_He3Pol", &T_err_He3Pol, "err_He3Pol/D");
  //double T_beamPol;	 Tout->Branch("beamPol", &T_beamPol, "beamPol/D");
  // Do the energy loss calculation here ...........

  // HCAL cut definitions
  double sbs_kick = kin_info.sbs_kick;
  vector<double> dx_p = kin_info.dx_p;
  vector<double> dy_p = kin_info.dy_p;
  double Nsigma_cut_dx_p = kin_info.Nsigma_dx_p;
  double Nsigma_cut_dy_p = kin_info.Nsigma_dy_p;
  vector<double> dx_n = kin_info.dx_n;
  vector<double> dy_n = kin_info.dy_n;
  double Nsigma_cut_dx_n = kin_info.Nsigma_dx_n;
  double Nsigma_cut_dy_n = kin_info.Nsigma_dy_n;
  vector<double> coin_time_cut = kin_info.coin_time_cut;
  double Nsigma_coin_time = kin_info.Nsigma_coin_time;
  vector<double> hcal_active_area = cut::hcal_active_area_data(); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p[1], dx_n[1], dy_p[1], hcal_active_area);

  // elastic cut limits
  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = kin_info.hcal_voffset;
  double hcal_hoffset = kin_info.hcal_hoffset;
  vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
  long nevent = 0, nevents = C->GetEntries(); 
  int treenum = 0, currenttreenum = 0, currentrunnum = 0;
  int IHWP_run = -100;  
  time_t run_time_unix;  

  cout<<"Processing "<<nevents<<" events"<<endl;
  
  while (C->GetEntry(nevent++)) {

    // print progress 
    if( nevent % 1000 == 0 ){ 
	//std::cout <<"theta_pq : "<<theta_pq<<endl;
	std::cout << nevent*100.0/nevents << "% \r";
    	std::cout.flush();
    }

    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();

      //Get the run number
      string s = C->GetFile()->GetName();
      int start = s.find("_stream0");
      start -= 4;
      int end = start + 4;
      T_runnum = stoi(s.substr(start,end - start));

      //Get the time this run started
      auto* Run_Data = C->GetFile()->Get<THaRunBase>("Run_Data");
      TDatime run_time = Run_Data->GetDate();
      run_time.Set(run_time.GetYear(),run_time.GetMonth(),run_time.GetDay(),run_time.GetHour(),run_time.GetMinute(),0);
      run_time_unix = run_time.Convert();

      // We must loop over a small subset to get the correct IHWP state for the entire run
      //First we do some stuff to max sure we dont reach the file limit
      int file_nevents = C->GetTree()->GetEntries();
      int max_events = 8000;
      if(file_nevents < max_events) max_events = file_nevents - 100;
      
      int start_event = nevent;
      while (C->GetEntry(start_event++) && start_event < nevent + max_events){
	if(IHWP == 1 || IHWP == -1) IHWP_run = IHWP;
      }
      C->GetEntry(nevent - 1);

    }
    
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;  
    
    if (!passedgCut) continue;
   
    //cout<<"global passed "<<endl;
    
    double bbcal_trig_time=0., hcal_trig_time=0.;
    for(int ihit=0; ihit<tdcElemN; ihit++){
      if(tdcElem[ihit]==5) bbcal_trig_time=tdcTrig[ihit];
      if(tdcElem[ihit]==0) hcal_trig_time=tdcTrig[ihit];
    }
 
    //Calculate absolute time of this event
    double time_interval = 4;  //in ns
    int time_rel = evtime*time_interval*1e-9 / 60; //in min, rounded

    TDatime time_abs(run_time_unix + time_rel * 60); //Add the relative minutes

    auto it = DBInfo.He3Pol.find(time_abs); // Get Polarization from the table
    if(it == DBInfo.He3Pol.end()){
      T_He3Pol = -1;
      T_err_He3Pol = -1;
    }
    else{
      T_He3Pol = it->second.first;
      T_err_He3Pol = it->second.second;
    }
    T_datetime = time_abs;
    
    //Timing Information
    T_hcal_time = atimeHCAL[0];
    T_bbcal_time = atimeSH;
      
    T_nhodo_clus = nhodo_clus;
    for(int iclus = 0; iclus < nhodo_clus; iclus++)
      T_hodo_time[iclus] = hodo_time[iclus];
      
    double coin_time = atimeHCAL[0] - atimeSH;  
    T_coin_time = coin_time; 
    h_coin_time->Fill(coin_time);
      
    coinCut = coin_time > coin_time_cut[0] - Nsigma_coin_time*coin_time_cut[1] && coin_time < coin_time_cut[0] + Nsigma_coin_time*coin_time_cut[1];
      
    //Grinch info
    T_grinch_track = grinch_track;
    T_grinch_clus_size = grinch_clus_size;
    
    //trigbits
    T_trigbits = trigbits;

    //Beam helicity information
    T_helicity = helicity;
    T_IHWP = IHWP_run;

    //BPMs and Rasters
    T_BPMAx = BPMAx;
    T_BPMAy = BPMAy;
    T_rawcurx = rawcurx;
    T_rawcury = rawcury;
    T_rawcur2x = rawcur2x;
    T_rawcur2y = rawcur2y;
    
    // kinematic parameters
    double ebeam = sbsconf.GetEbeam();       // Expected beam energy (GeV) [Get it from EPICS, eventually]
    double ebeam_corr = ebeam; //- MeanEloss;
    double precon = p[0]; //+ MeanEloss_outgoing

    // constructing the 4 vectors
    // Reaction    : e + e' -> N + N'
    // Conservation: Pe + Peprime = PN + PNprime 
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0,0,ebeam_corr,ebeam_corr);   // incoming e-
    TLorentzVector Peprime(px[0] * (precon/p[0]),   // scattered e-
			   py[0] * (precon/p[0]),
			   pz[0] * (precon/p[0]),
			   precon);                 
    TLorentzVector PN;                              // target nucleon [Ntype ??]
    kine::SetPN(Ntype, PN);

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pcentral = kine::pcentral(ebeam_corr, etheta, Ntype.Data());

    double nu = 0.;                   // energy of the virtual photon
    double pN_expect = 0.;            // expected recoil nucleon momentum
    double thetaN_expect = 0.;        // expected recoil nucleon theta
    double phiN_expect = ephi + constant::pi; 
    // Different modes of calculation. Goal is to achieve the best resolution
    // model 0 = uses reconstructed p as independent variable
    // model 1 = uses reconstructed angles as independent variable 
    // model 2 = uses 4-vector calculation 
    TVector3 pNhat;                   // 3-momentum of the recoil nucleon (Unit)
    double Q2recon, W2recon;
    if (model == 0) {
      nu = Pe.E() - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype.Data());
      thetaN_expect = acos((Pe.E() - Peprime.Pz()) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());
    } else if (model == 1) {
      nu = Pe.E() - pcentral;
      pN_expect = kine::pN_expect(nu, Ntype.Data());
      thetaN_expect = acos((Pe.E() - pcentral*cos(etheta)) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());
    } else if (model == 2) {
      TLorentzVector q = Pe - Peprime; // 4-momentum of virtual photon
      nu = q.E();
      pNhat = q.Vect().Unit();
      pN_expect = q.P();//spatial magnitude
      Q2recon = -q.M2();
      W2recon = (PN + q).M2();
      thetaN_expect = q.Theta();
    }
    h_Q2->Fill(Q2recon); 
    double Wrecon = sqrt(max(0., W2recon));
    double dpel = Peprime.E()/pcentral - 1.0; 
    h_dpel->Fill(dpel);

    
    TVector3 enhat_tgt( thtgt[0], phtgt[0], 1.0 );
    enhat_tgt = enhat_tgt.Unit();
	
    TVector3 enhat_fp( thfp[0], phfp[0], 1.0 );
    enhat_fp = enhat_fp.Unit();
	
    TVector3 enhat_fp_rot = enhat_fp.X() * GEMxaxis + enhat_fp.Y() * GEMyaxis + enhat_fp.Z() * GEMzaxis;

    double thetabend = acos( enhat_fp_rot.Dot( enhat_tgt ) );
    
    double pphi = atan(py_sbs[0]/px_sbs[0]) ;
    double ptheta = acos(pz_sbs[0]/p_sbs[0]);

    T_ebeam = Pe.E();

    //cout<<"here 0"<<endl;

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pcentral = pcentral;
    T_pN_expect = pN_expect;
    T_ptheta_cal = thetaN_expect;
    T_pphi_cal = phiN_expect;
    T_ptheta = ptheta;
    T_pphi = pphi;
    
    T_ntrack = ntrack;
    T_ntrack_hits = bb_gem_nhits[0];
    T_vz = vz[0];
    T_vx = vx[0];
    T_vy = vy[0];
    T_xtgt = xtgt[0];
    T_ytgt = ytgt[0];
    T_thtgt = thtgt[0];
    T_phtgt = phtgt[0];
    T_thetabend = thetabend;
    T_xfp = xfp[0];
    T_yfp = yfp[0];
    T_thfp = thfp[0];
    T_phfp = phfp[0];
    T_trP = p[0];
    T_trPx = px[0];
    T_trPy = py[0];
    T_trPz = pz[0];
    
    T_ntrack_sbs = ntrack_sbs; 
    T_ntrack_hits_sbs = sbs_gem_nhits[0];
    T_vz_sbs = vz_sbs[0];
    T_vx_sbs = vx_sbs[0];
    T_vy_sbs = vy_sbs[0];
    //T_xtgt_sbs = xtgt_sbs[0];
    T_ytgt_sbs = ytgt_sbs[0];
    T_thtgt_sbs = thtgt_sbs[0];
    T_phtgt_sbs = phtgt_sbs[0];
    //T_thetabend_sbs = thetabend_sbs;
    T_xfp_sbs = xfp_sbs[0];
    T_yfp_sbs = yfp_sbs[0];
    T_thfp_sbs = thfp_sbs[0];
    T_phfp_sbs = phfp_sbs[0];
    T_trP_sbs = p_sbs[0];
    T_trPx_sbs = px_sbs[0];
    T_trPy_sbs = py_sbs[0];
    T_trPz_sbs = pz_sbs[0];
    T_trx_sbs = xTr_sbs[0];
    T_try_sbs = yTr_sbs[0];
    T_trth_sbs = thTr_sbs[0];
    T_trph_sbs = phTr_sbs[0];

    T_ePS = ePS;
    T_xPS = xPS;
    T_eSH = eSH;
    T_xSH = xSH;
    T_ySH = ySH;

    T_eHCAL = eHCAL[0];
    //T_e_cHCAL = e_cHCAL[0];
    T_xHCAL = xHCAL[0];
    T_yHCAL = yHCAL[0];

    T_nclus_HCAL = nclusHCAL;
    T_nblk_HCAL = nblkHCAL;

    //cout<<"here 00"<<endl;
    /*
    for (int i = 0; i<nclusHCAL; ++i){
      T_hcal_clus_E[i] = clusEHCAL[i];
      //T_hcal_clus_atime[i] = clusatimeHCAL[i];
      //T_hcal_clus_tdctime[i] = clustdctimeHCAL[i];
      T_hcal_clus_id[i] = blkid[i];
      T_hcal_clus_nblk[i] = nblk[i]; 
    }
    */

    //make the outfile too big
    /*for (int i=0; i<nblkHCAL;++i){
      T_hcal_clus_mem_E[i] = clusmemEHCAL[i];
      //T_hcal_clus_mem_atime[i] = clusmematimeHCAL[i];
      //T_hcal_clus_mem_tdctime[i] = clusmemtdctimeHCAL[i];
      T_hcal_clus_mem_id[i] = clusmemblkid[i];
    }
    */

    //cout<<"here 3"<<endl;

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    double theta_pq = kine::theta_pq(vertex, HCAL_origin, pNhat, xHCAL[0], yHCAL[0]);
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    double dx = xHCAL[0] - xyHCAL_exp[0];  
    double dy = yHCAL[0] - xyHCAL_exp[1]; 

    //cout<<"here 2"<<endl;
    //if( nevent % 1000 == 0 ){
       // std::cout <<"HCAL x,y,z: "<<HCAL_origin[0]<<" "<<HCAL_origin[1]<<" "<<HCAL_origin[2]<<endl;
        //std::cout << nevent*100.0/nevents << "% \r";
        //std::cout.flush();
    //}

    T_theta_pq = theta_pq;
    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_dx = dx;
    T_dy = dy;

    // HCAL active area and safety margin cuts [Fiducial region]
    bool AR_cut = cut::inHCAL_activeA(xHCAL[0], yHCAL[0], hcal_active_area);
    bool FR_cut = cut::inHCAL_fiducial(xyHCAL_exp[0], xyHCAL_exp[1], sbs_kick, hcal_safety_margin);

    // HCAL cuts
    pCut = pow((dx-dx_p[0]) / (dx_p[1]*Nsigma_cut_dx_p), 2) + pow((dy-dy_p[0]) / (dy_p[1]*Nsigma_cut_dy_p), 2) <= 1.;
    nCut = pow((dx-dx_n[0]) / (dx_n[1]*Nsigma_cut_dx_n), 2) + pow((dy-dy_n[0]) / (dy_n[1]*Nsigma_cut_dy_n), 2) <= 1.;

    fiduCut = AR_cut && FR_cut;
    WCut = W2recon > W2min && W2recon < W2max;

    //cout<<"here 000"<<endl;

    // W cut and coincidence cut
    if (WCut && coinCut) {
      // fiducial cut
      //if (fiduCut) { //No fiducial cut for GEN?
	h_dxHCAL->Fill(dx);
	h_dyHCAL->Fill(dy);
	h2_rcHCAL->Fill(cblkHCAL[0], rblkHCAL[0]);
	h2_dxdyHCAL->Fill(dy, dx);
	
	if (pCut) h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
	if (nCut) h2_xyHCAL_n->Fill(xyHCAL_exp[1], xyHCAL_exp[0]);
	//}
    }
    h_W->Fill(Wrecon);
    // fiducial cut but no W cut
    //if (fiduCut) { //No fiducial cut for GEN
    //h_W_cut->Fill(Wrecon);
      if (pCut || nCut) { 
	h_W_cut->Fill(Wrecon);
      } else {
	h_W_acut->Fill(Wrecon);
      }
      //}
      
      //cout<<"here 1";

      Tout->Fill();
  } // event loop
  std::cout << std::endl << std::endl;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
  c1->Divide(2,2);

  c1->cd(1); h2_dxdyHCAL->Draw("colz");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0); Ep_p.SetLineColor(2); Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p[0], dx_p[0], Nsigma_cut_dy_p*dy_p[1], Nsigma_cut_dx_p*dx_p[1], 0,360,0);
  TEllipse Ep_n;
  Ep_n.SetFillStyle(0); Ep_n.SetLineColor(3); Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(dy_n[0], dx_n[0], Nsigma_cut_dy_n*dy_n[1], Nsigma_cut_dx_n*dx_n[1], 0,360,0);
 
  c1->cd(2);
  h_W->Draw(); h_W->SetLineColor(1);
  h_W_cut->Draw("same"); h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");

  c1->cd(3);
  h2_xyHCAL_p->Draw("colz");
  Utilities::DrawArea(hcal_active_area);
  c1->Write();
  
  fout->Write();
  sw->Delete();
   
  return 0;
}
