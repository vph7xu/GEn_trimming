//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   This script will "replay" the optics of a run so that we
//   do not need do a long full replay when testing a new 
//   optics calibration. It will take the analyzed root file 
//   from Quasielastic_ana.C and will copy all trees except
//   the optics. The optics will be calculated from the new
//   matrix which must be an input below in the script.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"
//#include "../../src/JSONManager.cpp"

void replay_run(const std::string cfg = "GEN3"){
  

  TFile *file;
  TFile *file_new;

  file = new TFile(Form("../outfiles/QE_data_%s_sbs100p_nucleon_p_model1_sbstrackingon.root",cfg.c_str()),"read");
  file_new = new TFile(Form("../outfiles/QE_new_matrix_%s_sbs100p_nucleon_p_model1_sbstrackingon.root",cfg.c_str()),"recreate");

   
  const std::string jmgr_file = "../../config/" + cfg + "_H2.cfg";
  if(cfg == "GEN2") const std::string jmgr_file = "../../config/" + cfg + "_H2_SBS100.cfg";
  

  JSONManager *jmgr = new JSONManager(jmgr_file);

  // seting up the desired SBS configuration
  TString conf = jmgr->GetValueFromKey_str("GEN_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  
  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromKey_str("Ntype");
  
  double CentAngle = sbsconf.GetBBtheta_rad();

  double GEMpitch = 10.0*TMath::DegToRad();
  TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
  TVector3 GEMyaxis(0,1,0);
  TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();

  //Matrix file from new optics
  TString Matrix_filename = "optics_gen2.txt";
  //Coefficients from momentum fit
  double A = 0.2775655;
  double B = 0.974132;
  double A_vy = 0;
  double B_vy = 0;

  if( cfg == "GEN2"){
    
    Matrix_filename = "optics_gen2.txt";
  
    /*
    A = 0.26182725;
    B = 0.974132;
    */
      
    A = 0.263077985;
    B = 0.999237626;
    A_vy = -0.0147632089;
    B_vy =-15.893095;
    
    
  }
  if( cfg == "GEN3"){
    Matrix_filename = "optics_gen3.txt";
    
    A = 0.27765103;
    B = 0.932092801;
    A_vy = 0.0175826672;
    B_vy = -33.8073321;
    
  }
  if( cfg == "GEN4b"){
    //Matrix_filename = "optics_gen4b.txt";
    Matrix_filename = "optics_gen4b_new.txt";

    A = 0.28588162;
    B = 0.900212499;
        
  }

  ifstream Matrix_file(Matrix_filename);


  int row_M = 0, col_M = 9;
  Matrix_file >> row_M;
  TMatrixD M(row_M,col_M);
   
  //Read in matrix file
  for(int row=0; row<row_M; row++)
    for(int col=0; col<col_M; col++) 
      Matrix_file >> M(row,col);
      
  
  
  TTree *T_old = (TTree*)file->Get("Tout");
   
  /////////// Variable Definitions /////////////////////////////////
  /////////////////////////////////////////////////////////////////
  int runnum;
  bool WCut, pCut, nCut, fiduCut, coinCut;
  std::vector<std::string> runvar = {"runnum","WCut","pCut","nCut","fiduCut","coinCut"};
  std::vector<void*> runvar_mem = {&runnum,&WCut, &pCut, &nCut, &fiduCut, &coinCut};
  setrootvar::setbranch(T_old,"",runvar,runvar_mem);

  double ebeam,nu,Q2,W2,dpel,ephi,etheta,pcentral;
  std::vector<std::string> kinevar = {"ebeam","nu","Q2","W2","dpel","ephi","etheta","pcentral"};
  std::vector<void*> kinevar_mem = {&ebeam,&nu,&Q2,&W2,&dpel,&ephi,&etheta,&pcentral};
  setrootvar::setbranch(T_old,"",kinevar,kinevar_mem);

  double vz,vx,vy,xtgt,ytgt,thtgt,phtgt,xfp,yfp,xpfp,ypfp,trP,trPx,trPy,trPz,thetabend;
  std::vector<std::string> trackvar = {"vz","vx","vy","xtgt","ytgt","thtgt","phtgt","thetabend","xfp","yfp","thfp","phfp","trP","trPx","trPy","trPz"};
  std::vector<void*> trackvar_mem = {&vz,&vx,&vy,&xtgt,&ytgt,&thtgt,&phtgt,&thetabend,&xfp,&yfp,&xpfp,&ypfp,&trP,&trPx,&trPy,&trPz};
  setrootvar::setbranch(T_old,"",trackvar,trackvar_mem);

  double ePS,xPS,eSH,xSH,ySH,eHCAL,xHCAL,yHCAL,xHCAL_exp,yHCAL_exp,dx,dy;
  std::vector<std::string> calvar = {"ePS","xPS","eSH","xSH","ySH","eHCAL","xHCAL","yHCAL","xHCAL_exp","yHCAL_exp","dx","dy"};
  std::vector<void*> calvar_mem = {&ePS,&xPS,&eSH,&xSH,&ySH,&eHCAL,&xHCAL,&yHCAL,&xHCAL_exp,&yHCAL_exp,&dx,&dy};
  setrootvar::setbranch(T_old,"",calvar,calvar_mem);

  int helicity,IHWP;
  double coinT_trig,hcal_time,bbcal_time,BPMAx,BPMAy,Rasterx,Rastery;
  std::vector<std::string> othervar = {"helicity","IHWP","coinT_trig","hcal_time","bbcal_time","BPMAx","BPMAy","Rasterx","Rastery"};
  std::vector<void*> othervar_mem = {&helicity,&IHWP,&coinT_trig,&hcal_time,&bbcal_time,&BPMAx,&BPMAy,&Rasterx,&Rastery};
  setrootvar::setbranch(T_old,"",othervar,othervar_mem); 

  const int maxHit = 1000;
  const int maxClus = 1000;
  int ngrinch_hits, nhodo_clus;
  double grinch_x[maxHit], grinch_y[maxHit],grinch_time[maxHit],hodo_time[maxClus];
  std::vector<std::string> grinchvar = {"xGRINCH","yGRINCH","grinch_time"};
  std::vector<void*> grinchvar_mem = {grinch_x,&grinch_y,&grinch_time};
  setrootvar::setbranch(T_old,"",grinchvar,grinchvar_mem); 
  setrootvar::setbranch(T_old,"","ngrinch_hits",&ngrinch_hits);
  setrootvar::setbranch(T_old,"","hodo_time",&hodo_time); 
  setrootvar::setbranch(T_old,"","nhodo_clus",&nhodo_clus); 

  ////////////////////////////////////////////////////////////////////////////

  file_new->cd();
  
  //TTree *T_new = T_old->CloneTree(0);
  TTree *T_new = new TTree("Tout","Tout");

  
  // Defining interesting ROOT tree branches 
  int T_runnum;         T_new->Branch("runnum", &T_runnum, "runnum/I");
  //cuts
  bool T_WCut;            T_new->Branch("WCut", &T_WCut, "WCut/B");
  bool T_pCut;            T_new->Branch("pCut", &T_pCut, "pCut/B");
  bool T_nCut;            T_new->Branch("nCut", &T_nCut, "nCut/B");
  bool T_fiduCut;         T_new->Branch("fiduCut", &T_fiduCut, "fiduCut/B");
  bool T_coinCut;         T_new->Branch("coinCut", &T_coinCut, "coinCut/B");
  //
  double T_ebeam;       T_new->Branch("ebeam", &T_ebeam, "ebeam/D");
  //kine
  double T_nu;          T_new->Branch("nu", &T_nu, "nu/D");
  double T_Q2;          T_new->Branch("Q2", &T_Q2, "Q2/D");
  double T_W2;          T_new->Branch("W2", &T_W2, "W2/D");
  double T_dpel;        T_new->Branch("dpel", &T_dpel, "dpel/D");
  double T_ephi;        T_new->Branch("ephi", &T_ephi, "ephi/D");
  double T_etheta;      T_new->Branch("etheta", &T_etheta, "etheta/D");
  double T_pcentral;    T_new->Branch("pcentral", &T_pcentral, "pcentral/D");
  //track
  double T_vz;          T_new->Branch("vz", &T_vz, "vz/D");
  double T_vx;          T_new->Branch("vx", &T_vx, "vx/D");
  double T_vy;          T_new->Branch("vy", &T_vy, "vy/D");
  double T_xtgt;        T_new->Branch("xtgt", &T_xtgt, "xtgt/D");
  double T_ytgt;        T_new->Branch("ytgt", &T_ytgt, "ytgt/D");
  double T_thtgt;       T_new->Branch("thtgt", &T_thtgt, "thtgt/D");
  double T_phtgt;       T_new->Branch("phtgt", &T_phtgt, "phtgt/D");
  double T_thetabend;   T_new->Branch("thetabend", &T_thetabend, "thetabend/D");
  double T_xfp;         T_new->Branch("xfp", &T_xfp, "xfp/D");
  double T_yfp;         T_new->Branch("yfp", &T_yfp, "yfp/D");
  double T_thfp;        T_new->Branch("thfp", &T_thfp, "thfp/D");
  double T_phfp;        T_new->Branch("phfp", &T_phfp, "phfp/D");  
  double T_trP;         T_new->Branch("trP", &T_trP, "trP/D");
  double T_trPx;         T_new->Branch("trPx", &T_trPx, "trPx/D");
  double T_trPy;         T_new->Branch("trPy", &T_trPy, "trPy/D");
  double T_trPz;         T_new->Branch("trPz", &T_trPz, "trPz/D");
  
  //BBCAL
  double T_ePS;         T_new->Branch("ePS", &T_ePS, "ePS/D"); 
  double T_xPS;         T_new->Branch("xPS", &T_xPS, "xPS/D"); 
  double T_eSH;         T_new->Branch("eSH", &T_eSH, "eSH/D"); 
  double T_xSH;         T_new->Branch("xSH", &T_xSH, "xSH/D"); 
  double T_ySH;         T_new->Branch("ySH", &T_ySH, "ySH/D"); 
  //HCAL
  double T_eHCAL;       T_new->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); 
  double T_xHCAL;       T_new->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); 
  double T_yHCAL;       T_new->Branch("yHCAL", &T_yHCAL, "yHCAL/D"); 
  double T_xHCAL_exp;   T_new->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D"); 
  double T_yHCAL_exp;   T_new->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D"); 
  double T_dx;          T_new->Branch("dx", &T_dx, "dx/D"); 
  double T_dy;          T_new->Branch("dy", &T_dy, "dy/D");
  //GRINCH
  int T_ngrinch_hits;           T_new->Branch("ngrinch_hits", &T_ngrinch_hits, "ngrinch_hits/I");
  double T_grinch_x[maxHit];    T_new->Branch("xGRINCH", &T_grinch_x, "xGRINCH[ngrinch_hits]/D");
  double T_grinch_y[maxHit];    T_new->Branch("yGRINCH", &T_grinch_y, "yGRINCH[ngrinch_hits]/D");
  
  //Timing Information
  double T_coinT_trig;           T_new->Branch("coinT_trig", &T_coinT_trig, "coinT_trig/D");
  double T_hcal_time;            T_new->Branch("hcal_time", &T_hcal_time, "hcal_time/D"); 
  double T_bbcal_time;           T_new->Branch("bbcal_time", &T_bbcal_time, "bbcal_time/D");
  double T_grinch_time[maxHit];  T_new->Branch("grinch_time", &T_grinch_time, "grinch_time[ngrinch_hits]/D");
  int T_nhodo_clus;              T_new->Branch("nhodo_clus", &T_nhodo_clus, "nhodo_clus/I");  
  double T_hodo_time[maxClus];   T_new->Branch("hodo_time", &T_hodo_time, "hodo_time[nhodo_clus]/D");

  //BPM and Raster information
  double T_BPMAx;     T_new->Branch("BPMAx", &T_BPMAx, "BPMAx/D");
  double T_BPMAy;     T_new->Branch("BPMAy", &T_BPMAy, "BPMAy/D");
  double T_rawcurx;  T_new->Branch("Rasterx", &T_rawcurx, "Raster_curx/D");
  double T_rawcury;  T_new->Branch("Rastery", &T_rawcury, "Raster_cury/D");

  //Helicity information
  int T_helicity;         T_new->Branch("helicity", &T_helicity, "helicity/I");
  int T_IHWP;         T_new->Branch("IHWP", &T_IHWP, "IHWP/I");

  int nevent = 0;
  int nevents = T_old->GetEntries();

  //Loop over all events
  while(T_old->GetEntry(nevent++)){

    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent*100.0/nevents << "% \r";
    std::cout.flush();
    
    double xptgt_new = 0;
    double yptgt_new = 0;
    double ytgt_new = 0;
    double xtgt_new = -vy;
    double pthetabend_new = 0;
    double vz_new = 0;
    
    //Calculate new optics variables from matrix
    for( int iter=0; iter<2; iter++ ){
       
      xptgt_new = 0.0;
      yptgt_new = 0.0;
      ytgt_new = 0.0;
      pthetabend_new = 0.0;
 
      for (int row=0; row<row_M; row++){
	xptgt_new += M(row,0)*pow(xfp,M(row,4))*pow(yfp,M(row,5))*pow(xpfp,M(row,6))*pow(ypfp,M(row,7))*pow(xtgt_new,M(row,8));
	yptgt_new += M(row,1)*pow(xfp,M(row,4))*pow(yfp,M(row,5))*pow(xpfp,M(row,6))*pow(ypfp,M(row,7))*pow(xtgt_new,M(row,8));
	ytgt_new += M(row,2)*pow(xfp,M(row,4))*pow(yfp,M(row,5))*pow(xpfp,M(row,6))*pow(ypfp,M(row,7))*pow(xtgt_new,M(row,8));
	pthetabend_new += M(row,3)*pow(xfp,M(row,4))*pow(yfp,M(row,5))*pow(xpfp,M(row,6))*pow(ypfp,M(row,7))*pow(xtgt_new,M(row,8));
	
      }
	  
      //beam left:
      vz_new = -ytgt_new / (sin(CentAngle) + cos(CentAngle)*yptgt_new);
       
      xtgt_new = -vy - vz_new * cos(CentAngle) * xptgt_new;
      
    }
    
    TVector3 enhat_tgt( xptgt_new, yptgt_new, 1.0 );
    enhat_tgt = enhat_tgt.Unit();
	
    TVector3 enhat_fp( xpfp, ypfp, 1.0 );
    enhat_fp = enhat_fp.Unit();
	
    TVector3 enhat_fp_rot = enhat_fp.X() * GEMxaxis + enhat_fp.Y() * GEMyaxis + enhat_fp.Z() * GEMzaxis;

    double thetabend_new = acos( enhat_fp_rot.Dot( enhat_tgt ) );
    
    //Calculate new momentum from coefficients
    double p_new = A * (1 + B * xptgt_new) / thetabend_new - (A_vy + B_vy * vy);

    double pz_new = p_new*sqrt( 1.0/(xptgt_new*xptgt_new+yptgt_new*yptgt_new+1.) );    
    double px_new = pz_new * xptgt_new;
    double py_new = pz_new * yptgt_new;

    TVector3 pvect_BB = TVector3(px_new, py_new, pz_new);

    px_new = +pvect_BB.Z()*sin(CentAngle)+pvect_BB.Y()*cos(CentAngle);
    py_new = -pvect_BB.X();
    pz_new = pvect_BB.Z()*cos(CentAngle)-pvect_BB.Y()*sin(CentAngle);
    
    TLorentzVector Peprime(px_new,   // scattered e-
			   py_new,
			   pz_new,
			   p_new);
    TLorentzVector Pe(0,0,ebeam,ebeam);
    double etheta_new = kine::etheta(Peprime);
    double pcentral_new = kine::pcentral(ebeam,etheta_new,Ntype);
    
    double nu = Pe.E() - pcentral_new;
    double pN_expect = kine::pN_expect(nu, Ntype);
    
    double Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta_new);
    double W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
    
    //Everything else should be copies from the tree
    T_runnum = runnum;
    //cuts
    T_WCut = WCut;
    T_pCut = pCut;
    T_nCut = nCut;
    T_fiduCut = fiduCut;
    T_coinCut = coinCut;
    //
    T_ebeam = ebeam;
    //kine
    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon;
    T_dpel = Peprime.E()/pcentral_new - 1.0; 
    T_ephi = ephi;
    T_etheta = etheta_new;
    T_pcentral = pcentral_new;
    //track
    T_vz = vz_new;
    T_vx = vx;
    T_vy = vy;
    T_xtgt = xtgt_new;
    T_ytgt = ytgt_new;
    T_thtgt = xptgt_new;
    T_phtgt = yptgt_new;
    T_thetabend = thetabend;
    T_xfp = xfp;
    T_yfp = yfp;
    T_thfp = xpfp;
    T_phfp = ypfp;
    T_trP = p_new;
    T_trPx = px_new;
    T_trPy = py_new;
    T_trPz = pz_new;
  
    //BBCAL
    T_ePS = ePS;
    T_xPS = xPS;
    T_eSH = eSH;
    T_xSH = xSH;
    T_ySH = ySH;
    //HCAL
    T_eHCAL = eHCAL;
    T_xHCAL = xHCAL;
    T_yHCAL = yHCAL;
    T_xHCAL_exp = xHCAL_exp;
    T_yHCAL_exp = yHCAL_exp;
    T_dx = dx;
    T_dy = dy;
    //GRINCH
    T_ngrinch_hits = ngrinch_hits;
    for(int ihit = 0; ihit < ngrinch_hits; ihit++){
      T_grinch_x[ihit] = grinch_x[ihit];
      T_grinch_y[ihit] = grinch_y[ihit];
      T_grinch_time[ihit] = grinch_time[ihit];
    }
  
    //Timing Information
    T_coinT_trig = coinT_trig;
    T_hcal_time = hcal_time;
    T_bbcal_time = bbcal_time;
    
    T_nhodo_clus = nhodo_clus;
    for(int ihit = 0; ihit < nhodo_clus; ihit++){
      T_hodo_time[ihit] = hodo_time[ihit];
    }

    //BPM and Raster information
    T_BPMAx = BPMAx;
    T_BPMAy = BPMAy;
    T_rawcurx = Rasterx;
    T_rawcury = Rastery;

    //Helicity information
    T_helicity = helicity;
    T_IHWP = IHWP;


    
    T_new->Fill();

  }

  file_new->cd();
  T_new->Write();
  file_new->Close();

}
