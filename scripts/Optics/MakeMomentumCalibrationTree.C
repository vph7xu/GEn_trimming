//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is find momemtum calibration
//   coefficients from H2 data. Input root file must be the
//   output from Quasielastic_ana.C
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//#include "GEM_cosmic_tracks.C"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TEventList.h"
#include "TCut.h"
#include <iostream>
#include <fstream>
#include "TMinuit.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TEventList.h"
#include "TLorentzVector.h" 
#include "TProfile.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TEllipse.h"
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include "TTreeFormula.h"

#include "../../include/gen-ana.h"
//#include "../../src/JSONManager.cpp"

double PI = TMath::Pi();

double Mp = 0.938272;
double Mn = 0.939565;
double clight = 2.99792458e-1; //m/ns

void draw_corr(TH2F *hist_orig, double (&fitresult)[2], TH2F *hist_new = NULL, vector<double> lincuts = {}){

  TH2F *hist = hist_orig;
  TH2F *hist_yrange = (TH2F*)hist->Clone("temp");
  
  
  double xmin = hist->GetXaxis()->GetXmin();
  double xmax = hist->GetXaxis()->GetXmax();
  
  hist_yrange->Reset();

  if(lincuts.size() != 0){

    double slope_low = (lincuts[1] - lincuts[0])/(xmax - xmin);
    double slope_high = (lincuts[3] - lincuts[2])/(xmax - xmin);
    double b_low = lincuts[1] - slope_low*xmax;
    double b_high = lincuts[3] - slope_high*xmax;

    for (Int_t binX = 1; binX <= hist->GetNbinsX(); binX++)
      {
	for (Int_t binY = 1; binY <= hist->GetNbinsY(); binY++)
	  {
	    double y = hist->GetYaxis()->GetBinCenter(binY);
	    double x = hist->GetXaxis()->GetBinCenter(binX);
	    if(y > slope_low*x + b_low && y < slope_high*x + b_high)
	      hist_yrange->SetBinContent(binX,binY,hist->GetBinContent(binX,binY));
	    
	  }
      }
  
    TLine *lin_low = new TLine(xmin,lincuts[0],xmax,lincuts[1]);
    lin_low->SetLineColor(kRed);  
    lin_low->Draw("same");

    TLine *lin_high = new TLine(xmin,lincuts[2],xmax,lincuts[3]);
    lin_high->SetLineColor(kRed);  
    lin_high->Draw("same");

    hist = (TH2F*)hist_yrange->Clone("temp2");
  }
  
  TF1 *linfit = new TF1("f1","pol1",xmin, xmax);

  hist->Fit("f1","qN");
  linfit->Draw("same");

  fitresult[0] = linfit->GetParameter(0);
  fitresult[1] = linfit->GetParameter(1);
  
  double mean_X = hist->GetMean(1);
  double mean_Y = hist->GetMean(2);
  double std_dev_X = hist->GetRMS(1);
  double std_dev_Y = hist->GetRMS(2);

  // Step 2: Calculate the covariance and variances
  double covariance = 0.0;
  double variance_X = 0.0;
  double variance_Y = 0.0;

  for (Int_t binX = 1; binX <= hist_orig->GetNbinsX(); binX++)
    {
      for (Int_t binY = 1; binY <= hist_orig->GetNbinsY(); binY++)
        {
	  double x = hist_orig->GetXaxis()->GetBinCenter(binX);
	  double y = hist_orig->GetYaxis()->GetBinCenter(binY);

	  double binContent = hist_orig->GetBinContent(binX, binY);

	  covariance += (x - mean_X) * (y - mean_Y) * binContent;
	  variance_X += std::pow(x - mean_X, 2) * binContent;
	  variance_Y += std::pow(y - mean_Y, 2) * binContent;
	 
	  if(hist_new != NULL){
	    double y_new = y - (linfit->GetParameter(1)*x + linfit->GetParameter(0));
	    int binY_new = hist_new->GetYaxis()->FindBin(y_new);
	    hist_new->SetBinContent(binX,binY_new,binContent);
	  }
	}
    }

}

//Main function
void MakeMomentumCalibrationTree( TString cfg = "GEN4", const char *outputfilename="NewMomentumFit.root"){

  
  gStyle->SetOptFit();

  //Reads from analyzed root file 
  TFile *Data_file;

  Data_file = new TFile("../outfiles/QE_new_matrix_" + cfg + "_sbs100p_nucleon_p_model1.root","read");

  TString jmgr_file = "../../config/" + cfg + "_H2.cfg";
  if(cfg == "GEN2") jmgr_file = "../../config/" + cfg + "_H2_SBS100.cfg";

  //JSONManager *jmgr = new JSONManager(jmgr_file);  

  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file,1);

  // seting up the desired SBS configuration
  //TString conf = jmgr->GetValueFromKey_str("GEN_config");
  //int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  //SBSconfig sbsconf(conf, sbsmag);

  TString conf = kin_info.conf;
  int sbsmag = kin_info.sbsmag;
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  //std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  
  
  //TString runinfo = "Runs: ";
  //TString runinfo2 = "";
  //for (int i=0; i<runnums.size(); i++) {
  //  if(i < 13) runinfo += Form("%i, ",runnums[i]);
  //  else runinfo2 += Form("%i, ",runnums[i]);
  //}

  //std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TCut globalcut = kin_info.globalcut;
  //TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);  
  
  //TString cutslist = "Cuts: " + gcut;

  //so as usual the main things to define are:
  // 1. List of files
  // 2. Global cuts
  // 3. Nominal beam energy (perhaps corrected for average energy loss along half target thickness
  // 4. File name for old momentum coefficients (optional?)
  // 5. BigBite central angle
  // 6. SBS central angle
  // 7. HCAL distance

  TTree *C = (TTree*)Data_file->Get("Tout");
  
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);  
  
  //Get all the coefficients from the SBS configuration
  int fitorder = 2;
  double ebeam=sbsconf.GetEbeam();
  double bbtheta = sbsconf.GetBBtheta_rad();
  double sbstheta = sbsconf.GetSBStheta_rad();
  double hcaldist = sbsconf.GetHCALdist(); //meters
  double BBdist = sbsconf.GetBBdist();
  double Ltgt = 0; //cm
  double rho_tgt = 0; //g/cc
  double rho_Al = 0; //g/cc

  double sbsdist = sbsconf.GetSBSdist();
  double sbsfield = sbsconf.GetSBSmag()*1.0/100; //fraction of max field:
  //double sbsmaxfield = 0.53/0.3; //About 1.77 Tesla.
  double sbsmaxfield = 1.3; //Tesla
  double Dgap = 48.0*2.54/100.0; //about 1.22 m
  
  double celldiameter = 1.6*2.54; //cm, right now this is a guess

  //Eventually we will grab HALLA:p from the EPICs tree for the beam energy:

  // 0.85 m deflection / 3.1 nucleon momentum @70% of max field
  // p * thetabend = 0.3 * BdL
  // deflection = tan(thetabend) * (11 - (2.25 + 0.6) )
    //BdL = p*thetabend/0.3 --> Bmax = p*thetabend/0.3/dL/0.7

    
  double Ztgt = 1.0;
  double Atgt = 1.0;
  double Mmol_tgt = 1.008; //g/mol
  
  //For energy-loss correction to beam energy:
  
  double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy:
  double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV:

  double uwallthick_LH2 = 0.0145; //cm
  double dwallthick_LH2 = 0.015; //cm
  double cellthick_LH2 = 0.02; //cm, this is a guess;

  int useAlshield = 0; //Use 1/8" aluminum shield on scattering chamber exit?

  double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch
  
  
  //Mean energy loss of the beam prior to the scattering:
  double MeanEloss = Ltgt/2.0 * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV:

  double MeanEloss_outgoing = celldiameter/2.0/sin(bbtheta) * rho_tgt * dEdx_tgt; //Approximately 1 MeV
  
  int usehcalcut = 1; //default to false:


  //vector<double> dx_p; //jmgr->GetVectorFromKey<double>("dx_p", dx_p);
  //vector<double> dy_p; //jmgr->GetVectorFromKey<double>("dy_p", dy_p);
  
  vector<double> dx_p = kin_info.dx_p;
  vector<double> dy_p = kin_info.dy_p;

  double nsigma = 2.0;

  double GEMpitch = 10.0*TMath::DegToRad();

  double dpelmin_fit = -0.03;
  double dpelmax_fit = 0.03;
  double W2min_fit = 0.5;
  double W2max_fit = 2.5; 

  //The following are the positions of the "first" row and column from HCAL database (top right block as viewed from upstream)
  double xoff_hcal = kin_info.hcal_voffset; //jmgr->GetValueFromKey<double>("hcal_voffset");
  double yoff_hcal = kin_info.hcal_hoffset; //jmgr->GetValueFromKey<double>("hcal_hoffset");
  
  //By default we fix the zero-order coefficient for pth and the first-order pth coefficient for xfp:
  
  int fix_pth0_flag = 1;
  int fix_pthx_flag = 1;

  double pth0 = 0.275;
  double pthx = 0.102;
  
  int order = 2; 
  //int order_ptheta = 1; //default to first-order expansion for momentum, and 2nd-order for vertices, angles. 
  
  double A_pth0 = 0.275;
  double B_pth0 = 0.61;
  double C_pth0 = -0.074;

  
  
  TString fname_oldcoeffs = "newcoeffs.dat"; //mainly for merging angle/vertex reconstruction coefficients with momentum coefficients:

  //Look at plots to decide the correct range for these cuts
  double thtgtmin_fit = -0.15;
  double thtgtmax_fit = 0.05;
  if(cfg == "GEN2") thtgtmax_fit = -0.01;

  //double Wmin_treefill = 0;
  //double Wmax_treefill = 

  int hcal_coordflag = 0; //0 = "new" HCAL coordinate system that matches transport.
                          //1 = "old" HCAL coordinate system that doesn't match transport
  
  
  
  //Note that both of these calculations neglect the Aluminum end windows and cell walls:
  
  MeanEloss = Ltgt/2.0 * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //hopefully the user has provided these quantities with correct units
  
  MeanEloss_outgoing = celldiameter/2.0/sin(bbtheta) * rho_tgt * dEdx_tgt +
    cellthick_LH2/sin(bbtheta) * rho_Al * dEdx_Al;
  if( useAlshield != 0 ) MeanEloss_outgoing += Alshieldthick * rho_Al * dEdx_Al;

  cout << "mean eloss in aluminum shield = " << Alshieldthick * rho_Al * dEdx_Al << endl;
  
  cout << "use Aluminum shield = " << useAlshield << endl;
						     
						     
  
  cout << "Mean E loss (beam electron) = " << MeanEloss << endl;
  cout << "BigBite theta = " << bbtheta * TMath::RadToDeg() << endl;
  cout << "Mean E loss (outgoing electron) = " << MeanEloss_outgoing << endl;
  cout << "Ebeam corrected (GeV) = " << ebeam - MeanEloss << endl;
 
					   
					   
  

  C->SetBranchStatus("*",0);

  // track var
  double p,px,py,pz;
  double vx,vy,vz;
  double xtgt,ytgt,thtgt,phtgt;
  double xfp,yfp,thfp,phfp;
  std::vector<std::string> trvar = {"trP","trPx","trPy","trPz","vx","vy","vz","xtgt","ytgt","thtgt","phtgt","xfp","yfp","thfp","phfp"};
  std::vector<void*> trvar_mem = {&p,&px,&py,&pz,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,&xfp,&yfp,&thfp,&phfp};
  setrootvar::setbranch(C,"",trvar,trvar_mem);
  

  // hcal clus var
  double EHCAL, xHCAL, yHCAL, dx, dy;
  std::vector<std::string> hcalclvar = {"eHCAL","xHCAL","yHCAL","dx","dy"};
  std::vector<void*> hcalclvar_mem = {&EHCAL,&xHCAL,&yHCAL,&dx,&dy};
  setrootvar::setbranch(C,"", hcalclvar, hcalclvar_mem);

  // bbcal sh var
  double ESH, xSH, ySH;
  std::vector<std::string> bbcalshclvar = {"eSH","xSH","ySH"}; 
  std::vector<void*> bbcalshclvar_mem = {&ESH,&xSH,&ySH}; 
  setrootvar::setbranch(C,"",bbcalshclvar,bbcalshclvar_mem);

  // bbcal ps var
  double EPS, xPS;
  std::vector<std::string> bbcalpsclvar = {"ePS","xPS"}; 
  std::vector<void*> bbcalpsclvar_mem = {&EPS,&xPS}; 
  setrootvar::setbranch(C,"",bbcalpsclvar,bbcalpsclvar_mem);

  // Kinematics
  double nu_recon, Q2recon, W2recon, dpel;
  std::vector<std::string> kinvar = {"nu","Q2","W2", "dpel"}; 
  std::vector<void*> kinvar_mem = {&nu_recon,&Q2recon,&W2recon,&dpel}; 
  setrootvar::setbranch(C,"",kinvar,kinvar_mem);

  // Cuts
  bool pCut, WCut;
  std::vector<std::string> cutvar = {"pCut","WCut"}; 
  std::vector<void*> cutvar_mem = {&pCut,&WCut}; 
  setrootvar::setbranch(C,"",cutvar,cutvar_mem);
  
  int nparams = 0;

  vector<int> xtar_expon;
  vector<int> xfp_expon;
  vector<int> xpfp_expon;
  vector<int> yfp_expon;
  vector<int> ypfp_expon;

  //ifstream foldcoeffs(fname_oldcoeffs.Data());
  
  for(int i=0; i<=order; i++){
    for(int j=0; j<=order-i; j++){
      for(int k=0; k<=order-i-j; k++){
  	for(int l=0; l<=order-i-j-k; l++){
  	  for(int m=0; m<=order-i-j-k-l; m++){
  	    nparams++;
  	    xtar_expon.push_back( i );
  	    xfp_expon.push_back( m );
  	    yfp_expon.push_back( l );
  	    xpfp_expon.push_back( k );
  	    ypfp_expon.push_back( j );
  	  }
  	}
      }
    }
  }

  //For our purposes here we only need to set up the augmented matrix for the p*thetabend matrix elements:
  TMatrixD M(nparams,nparams);
  TVectorD b_pth(nparams);

  for( int i=0; i<nparams; i++ ){
    for( int j=0; j<nparams; j++ ){
      M(i,j) = 0.0;
    }
    b_pth(i) = 0.0;
  }


  double pcentral = ebeam/(1.+ebeam/Mp*(1.-cos(bbtheta)));

  TFile *fout = new TFile(outputfilename,"RECREATE");

  //After we get this working, we'll also want correlation plots of dpel and W vs. focal plane and target variables, new and old:
  TH1D *hdpel_old = new TH1D("hdpel_old","Old MEs;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH1D *hW_old = new TH1D("hW_old","Old MEs;W (GeV);",250,0,2);
  TH1D *hp_old = new TH1D("hp_old","Old MEs;p (GeV);",250,0.5*pcentral,1.5*pcentral);
  TH1D *hpthetabend_old = new TH1D("hpthetabend_old","Old MEs; p#theta_{bend};",250,0.2,0.35);

  TH2D *hxpfp_vs_xfp_old = new TH2D("hxpfp_vs_xfp_old",";xfp (m); xpfp",250,-0.75,0.75,250,-0.4,0.4);
  
  TH2D *hdpel_xfp_old = new TH2D("hdpel_xfp_old", "Old MEs; xfp (m); p/p_{elastic}(#theta)-1",250,-0.75,0.75,250,-0.125,0.125);
  TH2D *hdpel_yfp_old = new TH2D("hdpel_yfp_old", "Old MEs; yfp (m); p/p_{elastic}(#theta)-1",250,-0.25,0.25,250,-0.125,0.125);
  TH2D *hdpel_xpfp_old = new TH2D("hdpel_xpfp_old", "Old MEs; xpfp; p/p_{elastic}(#theta)-1",250,-0.4,0.4,250,-0.125,0.125);
  TH2D *hdpel_ypfp_old = new TH2D("hdpel_ypfp_old", "Old MEs; ypfp; p/p_{elastic}(#theta)-1",250,-0.15,0.15,250,-0.125,0.125);

  TH2D *hpthetabend_xfp_old = new TH2D("hpthetabend_xfp_old", "Old MEs; xfp (m); p#theta_{bend}",250,-0.75,0.75,250,0.2,0.35);
  TH2D *hpthetabend_yfp_old = new TH2D("hpthetabend_yfp_old", "Old MEs; yfp (m); p#theta_{bend}",250,-0.25,0.25,250,0.2,0.35);
  TH2D *hpthetabend_xpfp_old = new TH2D("hpthetabend_xpfp_old", "Old MEs; xpfp; p#theta_{bend}",250,-0.4,0.4,250,0.2,0.35);
  TH2D *hpthetabend_ypfp_old = new TH2D("hpthetabend_ypfp_old", "Old MEs; ypfp; p#theta_{bend}",250,-0.15,0.15,250,0.2,0.35);

  TH2D *hpthetabend_xptar_old = new TH2D("hpthetabend_xptar_old", "Old MEs; xptar; p#theta_{bend}",250,-0.4,0.4,250,0.0,0.5);
  TH2D *hpthetabend_yptar_old = new TH2D("hpthetabend_yptar_old", "Old MEs; yptar; p#theta_{bend}",250,-0.15,0.15,250,0.2,0.35);
  TH2D *hpthetabend_ytar_old = new TH2D("hpthetabend_ytar_old", "Old MEs; ytar; p#theta_{bend}",250,-0.15,0.15,250,0.2,0.35);

  TH2D *hdpel_xptar_old = new TH2D("hdpel_xptar_old", "Old MEs; xptar; p/p_{elastic}(#theta)-1",250,-0.4,0.4,250,-0.125,0.125);
  TH2D *hdpel_yptar_old = new TH2D("hdpel_yptar_old", "Old MEs; yptar; p/p_{elastic}(#theta)-1",250,-0.15,0.15,250,-0.125,0.125);
  TH2D *hdpel_ytar_old = new TH2D("hdpel_ytar_old", "Old MEs; ytar; p/p_{elastic}(#theta)-1",250,-0.15,0.15,250,-0.125,0.125);

  TH1D *hdx_HCAL_old = new TH1D("hdx_HCAL_old","Old MEs; xHCAL - xBB (m);", 250, -2.5,2.5);
  TH1D *hdy_HCAL_old = new TH1D("hdy_HCAL_old","Old MEs; yHCAL - yBB (m);", 250, -1.25,1.25);
  
  TH2D *hdxdy_HCAL_old = new TH2D("hdxdy_HCAL_old","Old MEs; yHCAL - yBB (m); xHCAL - xBB (m)", 250, -3.25,3.25, 250,-3.5,3.5);
  TH2D *hdxdy_HCAL_new = new TH2D("hdxdy_HCAL_new","Old MEs; yHCAL - yBB (m); yHCAL - yBB (m)", 250, -3.25,3.25, 250,-3.5,3.5);
  
  TH1D *hW_new = new TH1D("hW_new","New MEs;W (GeV);",250,0,2);
  TH1D *hdpel_new = new TH1D("hdpel_new","New MEs;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH2F *hhp_vs_vy = new TH2F("hhp_vs_vy","Momentum Correlation;vy;p",100,-0.004,0.004,100,1,4);
  TH2F *hhdpel_vs_vy = new TH2F("hhdpel_vs_vy","Momentum Correlation;vy;p/p_{elastic}(#theta)-1",100,-0.004,0.004,100,-0.1,0.1);
  TH2F *hhdpel_vs_vy_temp = new TH2F("hhdpel_vs_vy_temp","Momentum Correlation;vy;p/p_{elastic}(#theta)-1",100,-0.004,0.004,100,-0.1,0.1);
  TH2F *hhdpel_vs_vy_new = new TH2F("hhdpel_vs_vy_new","Momentum Correlation New;vy;p/p_{elastic}(#theta)-1",100,-0.004,0.004,100,-0.1,0.1);
  TH2F *hhdpel_vs_xtar = new TH2F("hhdpel_vs_xtar", "Momentum Correlation;x_{tg};p/p_{elastic}(#theta)-1",50,-0.004,0.004,50,-0.1,0.1);
  TH2F *hhdpel_vs_xtar_temp = new TH2F("hhdpel_vs_xtar_temp", "Momentum Correlation;x_{tg};p/p_{elastic}(#theta)-1",50,-0.004,0.004,50,-0.1,0.1);
  TH2F *hhdpel_vs_xtar_new = new TH2F("hhdpel_vs_xtar_new", "Momentum Correlation New;x_{tg};p/p_{elastic}(#theta)-1",50,-0.004,0.004,50,-0.1,0.1);
  

  TH1D *hW_new_nocut = new TH1D("hW_new_nocut","New MEs;W (GeV);",250,0,2);
  TH1D *hdpel_new_nocut = new TH1D("hdpel_new_nocut","New MEs;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  
  TTree *Tout = new TTree("Tout","Tree containing variables for momentum calibration");

  double T_ebeam, T_etheta, T_ephi, T_precon, T_pelastic, T_thetabend, T_dpel, T_W2;
  double T_pincident;
  double T_xfp, T_yfp, T_thfp, T_phfp;
  double T_thtgt, T_phtgt, T_ytgt, T_xtgt;
  double T_vx, T_vy, T_vz;
  double T_BBdist, T_BBtheta;
  double T_HCALdist, T_HCALtheta;
  double T_xHCAL, T_yHCAL, T_EHCAL, T_deltax, T_deltay;
  double T_xHCAL_expect, T_yHCAL_expect;
  double T_pp_expect, T_ptheta_expect, T_pphi_expect;
  double T_EPS, T_ESH, T_Etot;
  double T_xSH, T_ySH;
  double T_Q2;
  double T_Q2_4vect, T_W2_4vect;
  double T_thetaq, T_phiq, T_px, T_py, T_pz, T_qx, T_qy, T_qz, T_qmag;
  double T_theta_recon_n, T_phi_recon_n, T_theta_recon_p, T_phi_recon_p; //Reconstructed nucleon angles under proton and neutron hypothesis
  double T_thetapq_n, T_phipq_n, T_thetapq_p, T_phipq_p;
  double T_dx_4vect, T_dy_4vect; //Here we want to use the 4-vector momentum transfer to calculate dx/dy
  double T_hodotime, T_hcaltime;
  double T_Lneutron, T_nTOFexpect;
  
  int HCALcut;
  int BBcut;
  int T_WCut;
  double ep_incident_avg = 0;
  int n_ep_inc = 0;
  
  Tout->Branch( "HCALcut", &HCALcut, "HCALcut/I");
  Tout->Branch( "BBcut", &BBcut, "BBcut/I");
  Tout->Branch( "WCut", &T_WCut, "WCut/I");
  Tout->Branch( "Ebeam", &T_ebeam, "Ebeam/D" );
  Tout->Branch( "Q2", &T_Q2, "Q2/D"); 
  Tout->Branch( "etheta", &T_etheta, "etheta/D");
  Tout->Branch( "ephi", &T_ephi, "ephi/D");
  Tout->Branch( "ep_recon", &T_precon, "ep_recon/D");
  Tout->Branch( "ep_elastic", &T_pelastic, "ep_elastic/D");
  Tout->Branch( "ep_incident", &T_pincident, "ep_incident/D");
  Tout->Branch( "thetabend", &T_thetabend, "thetabend/D");
  Tout->Branch( "dpel", &T_dpel, "dpel/D");
  Tout->Branch( "W2", &T_W2, "W2/D");
  Tout->Branch( "xfp", &T_xfp, "xfp/D");
  Tout->Branch( "yfp", &T_yfp, "yfp/D");
  Tout->Branch( "thfp", &T_thfp, "thfp/D");
  Tout->Branch( "phfp", &T_phfp, "phfp/D");
  Tout->Branch( "thtgt", &T_thtgt, "thtgt/D");
  Tout->Branch( "phtgt", &T_phtgt, "phtgt/D");
  Tout->Branch( "ytgt", &T_ytgt, "ytgt/D");
  Tout->Branch( "xtgt", &T_xtgt, "xtgt/D");
  Tout->Branch( "vx", &T_vx, "vx/D");
  Tout->Branch( "vy", &T_vy, "vy/D");
  Tout->Branch( "vz", &T_vz, "vz/D");
  Tout->Branch( "BBdist", &T_BBdist, "BBdist/D");
  Tout->Branch( "BBtheta", &T_BBtheta, "BBtheta/D");
  Tout->Branch( "HCALdist", &T_HCALdist, "HCALdist/D");
  Tout->Branch( "HCALtheta", &T_HCALtheta, "HCALtheta/D");
  Tout->Branch( "xHCAL", &T_xHCAL, "xHCAL/D");
  Tout->Branch( "yHCAL", &T_yHCAL, "yHCAL/D");
  Tout->Branch( "xHCAL_expect", &T_xHCAL_expect, "xHCAL_expect/D");
  Tout->Branch( "yHCAL_expect", &T_yHCAL_expect, "yHCAL_expect/D");
  Tout->Branch( "EHCAL", &T_EHCAL, "EHCAL/D");
  Tout->Branch( "deltax", &T_deltax, "deltax/D");
  Tout->Branch( "deltay", &T_deltay, "deltay/D");
  Tout->Branch( "pp_expect", &T_pp_expect, "pp_expect/D");
  Tout->Branch( "ptheta_expect", &T_ptheta_expect, "ptheta_expect/D");
  Tout->Branch( "pphi_expect", &T_pphi_expect, "pphi_expect/D");
  Tout->Branch( "EPS", &T_EPS, "EPS/D");
  Tout->Branch( "ESH", &T_ESH, "ESH/D");
  Tout->Branch( "Etot", &T_Etot, "Etot/D");
  Tout->Branch( "xSH", &T_xSH, "xSH/D");
  Tout->Branch( "ySH", &T_ySH, "ySH/D");
  Tout->Branch( "Q2_4vect", &T_Q2_4vect, "Q2_4vect/D");
  Tout->Branch( "W2_4vect", &T_W2_4vect, "W2_4vect/D");
  Tout->Branch( "epx", &T_px, "epx/D");
  Tout->Branch( "epy", &T_py, "epy/D");
  Tout->Branch( "epz", &T_pz, "epz/D");
  Tout->Branch( "qx", &T_qx, "qx/D" );
  Tout->Branch( "qy", &T_qy, "qy/D" );
  Tout->Branch( "qz", &T_qz, "qz/D" );
  Tout->Branch( "qmag", &T_qmag, "qmag/D" );
  Tout->Branch( "thetaq", &T_thetaq, "thetaq/D");
  Tout->Branch( "phiq", &T_phiq, "phiq/D" );
  Tout->Branch( "thetarecon_n", &T_theta_recon_n, "thetarecon_n/D" );
  Tout->Branch( "phirecon_n", &T_phi_recon_n, "phirecon_n/D" );
  Tout->Branch( "thetarecon_p", &T_theta_recon_p, "thetarecon_p/D" );
  Tout->Branch( "phirecon_p", &T_phi_recon_p, "phirecon_p/D" );
  Tout->Branch( "thetapq_n", &T_thetapq_n, "thetapq_n/D" );
  Tout->Branch( "thetapq_p", &T_thetapq_p, "thetapp_n/D" );
  //Tout->Branch( "phipq_n", &T_phipq_n, "phipq_n/D" );
  //Tout->Branch( "phipq_p", &T_phipq_p, "phipp_n/D" );
  Tout->Branch( "deltax_4vect", &T_dx_4vect, "deltax_4vect/D" );
  Tout->Branch( "deltay_4vect", &T_dy_4vect, "deltay_4vect/D" );
  Tout->Branch( "hodotime", &T_hodotime, "hodotime/D");
  Tout->Branch( "hcaltime", &T_hcaltime, "hcaltime/D");
  Tout->Branch( "Lneutron", &T_Lneutron, "Lneutron/D");
  Tout->Branch( "nTOFexpect", &T_nTOFexpect, "nTOFexpect/D");

  long ntotal = C->GetEntries();
  
  long nevent=0;

  //First pass: accumulate sums required for the fit: 
  
  int treenum=0, currenttreenum=0; 

  while( C->GetEntry(nevent++)){
    if( nevent % 100000 == 0 ) std::cout << nevent*100.0/ntotal << "% \r";
    std::cout.flush();      
    


    
      //The first thing we want to do is to calculate the "true" electron momentum incident on BigBite:
      

      double Ebeam_corrected = ebeam - MeanEloss;

      T_ebeam = Ebeam_corrected;
      
      double etheta = acos(pz/p);

      T_etheta = etheta;
      
      // Calculate the expected momentum of an elastically scattered electron at the reconstructed scattering angle and then correct it for the mean energy loss of the
      // electron on its way out of the target:
      double pelastic = Ebeam_corrected/(1.+Ebeam_corrected/Mp*(1.0-cos(etheta))); 
      
      double precon = p + MeanEloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)

      T_precon = precon;
      T_pelastic = pelastic; 
      T_pincident = pelastic - MeanEloss_outgoing;

      T_Q2 = Q2recon;
      T_W2 = W2recon;

      T_xfp = xfp;
      T_yfp = yfp;
      T_thfp = thfp;
      T_phfp = phfp;

      T_thtgt = thtgt;
      T_phtgt = phtgt;
      T_ytgt = ytgt;
      T_xtgt = xtgt;

      T_vx = vx;
      T_vy = vy;
      T_vz = vz;

      T_BBdist = BBdist;
      T_BBtheta = bbtheta;
      T_HCALdist = hcaldist;
      T_HCALtheta = sbstheta;

      
      double Wrecon = sqrt( std::max(0.0,W2recon) );

      T_dpel = dpel;
      
      bool passed_HCAL_cut = true;

      TVector3 vertex( 0, 0, vz );
      TLorentzVector Pbeam(0,0,Ebeam_corrected,Ebeam_corrected);
      TLorentzVector Kprime(precon*px/p,precon*py/p,precon*pz/p,precon);
      TLorentzVector Ptarg(0,0,0,Mp);
      
      TLorentzVector q = Pbeam - Kprime;

      T_Q2_4vect = -q.M2();
      T_W2_4vect = (Ptarg + q).M2();

      TVector3 qvect = q.Vect();
      T_thetaq = qvect.Theta();
      T_phiq = qvect.Phi();
      T_px = Kprime.Px();
      T_py = Kprime.Py();
      T_pz = Kprime.Pz();

      T_qx = q.Px();
      T_qy = q.Py();
      T_qz = q.Pz();

      T_qmag = qvect.Mag();
      
      double ephi = atan2( py, px );
      //double etheta = acos( pz/p );

      T_ephi = ephi;
      
      //      double nu = Ebeam_corrected - p;

      // Usually when we are running this code, the angle reconstruction is already well calibrated, but the momentum reconstruction is
      // unreliable; use pel(theta) as electron momentum for kinematic correlation:
      double nu = Ebeam_corrected - pelastic; 
      
      double pp_expect = sqrt(pow(nu,2)+2.*Mp*nu); //sqrt(nu^2 + Q^2) = sqrt(nu^2 + 2Mnu)
      double pphi_expect = ephi + PI;
      //double ptheta_expect = acos( (Ebeam_corrected-pz)/pp_expect ); //will this give better resolution than methods based on electron angle only? Not entirely clear

      T_pp_expect = pp_expect;
      T_pphi_expect = pphi_expect;
      
      double ptheta_expect = acos( (Ebeam_corrected-pelastic*cos(etheta))/pp_expect ); //will this give better resolution than methods based on electron angle only? Not entirely clear

      T_ptheta_expect = ptheta_expect;

      T_EPS = EPS;
      T_ESH = ESH;
      T_Etot = EPS + ESH;
      T_xSH = xSH;
      T_ySH = ySH;
      
      TVector3 pNhat( sin(ptheta_expect)*cos(pphi_expect), sin(ptheta_expect)*sin(pphi_expect), cos(ptheta_expect) );

      TVector3 HCAL_zaxis(-sin(sbstheta),0,cos(sbstheta));
      TVector3 HCAL_xaxis(0,1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
      
      TVector3 HCAL_origin = hcaldist * HCAL_zaxis + xoff_hcal*HCAL_xaxis + yoff_hcal*HCAL_yaxis;

      
      double sintersect = (HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot(HCAL_zaxis));
      //Straight-line projection to the surface of HCAL:
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;

      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

      TVector3 qunit = qvect.Unit();

      double s4vect = (HCAL_origin - vertex).Dot( HCAL_zaxis )/ (qunit.Dot( HCAL_zaxis ) );
      TVector3 HCAL_intersect4 = vertex + s4vect * qunit;

      T_Lneutron = (HCAL_intersect4 - vertex).Mag();

      double betan_expect = T_qmag / sqrt( pow( T_qmag,2) + pow( Mn,2 ) );
      
      T_nTOFexpect = T_Lneutron/(betan_expect * clight); //should be in ns

      T_dx_4vect = xHCAL - (HCAL_intersect4 - HCAL_origin).Dot( HCAL_xaxis );
      T_dy_4vect = yHCAL - (HCAL_intersect4 - HCAL_origin).Dot( HCAL_yaxis );

      TVector3 HCALpos = HCAL_origin + xHCAL * HCAL_xaxis + yHCAL * HCAL_yaxis;

      //Calculate expected proton deflection using crude model:
      double BdL = sbsfield * sbsmaxfield * Dgap;

      //thetabend = 0.3 * BdL: 
      double proton_deflection = tan( 0.3 * BdL / qvect.Mag() ) * (hcaldist - (sbsdist + Dgap/2.0) );

      //cout << "Expected proton deflection = " << proton_deflection << " meters" << endl;
      
      TVector3 NeutronDirection = (HCALpos - vertex).Unit();
      T_theta_recon_n = acos( NeutronDirection.Z() );
      T_phi_recon_n = TMath::ATan2( NeutronDirection.Y(), NeutronDirection.X() );

      TVector3 ProtonDirection = (HCALpos + proton_deflection * HCAL_xaxis - vertex).Unit();
      T_theta_recon_p = acos( ProtonDirection.Z() );
      T_phi_recon_p = TMath::ATan2( ProtonDirection.Y(), ProtonDirection.X() );

      T_thetapq_n = acos( NeutronDirection.Dot( qvect.Unit() ) );
      T_thetapq_p = acos( ProtonDirection.Dot( qvect.Unit() ) );
      
      //Now we need to calculate the "true" trajectory bend angle for the electron from the reconstructed angles:
      TVector3 enhat_tgt( thtgt, phtgt, 1.0 );
      enhat_tgt = enhat_tgt.Unit();
	
      TVector3 enhat_fp( thfp, phfp, 1.0 );
      enhat_fp = enhat_fp.Unit();

      TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
      TVector3 GEMyaxis(0,1,0);
      TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();
	
      TVector3 enhat_fp_rot = enhat_fp.X() * GEMxaxis + enhat_fp.Y() * GEMyaxis + enhat_fp.Z() * GEMzaxis;

      double thetabend = acos( enhat_fp_rot.Dot( enhat_tgt ) );

      T_thetabend = thetabend;
      
      if( W2recon >= W2min_fit && W2recon <= W2max_fit ){
	hdx_HCAL_old->Fill( dx );
	hdy_HCAL_old->Fill( dy );
	hdxdy_HCAL_old->Fill( dy,
			      dx );
      }

      T_xHCAL = xHCAL;
      T_yHCAL = yHCAL;
      T_EHCAL = EHCAL;
      T_deltax = dx;
      T_deltay = dy;

      T_xHCAL_expect = xexpect_HCAL;
      T_yHCAL_expect = yexpect_HCAL;
      
      if( usehcalcut != 0 ){
	passed_HCAL_cut = pCut;
      }

      T_WCut = WCut;
      HCALcut = pow( (dx - dx_p[0])/dx_p[1], 2 ) +
	pow( (dy - dy_p[0])/dy_p[1], 2 ) <= pow(nsigma,2); 

      BBcut = W2recon >= W2min_fit && W2recon <= W2max_fit && dpel >= dpelmin_fit && dpel <= dpelmax_fit;
      
      if( passed_HCAL_cut ){

	hdpel_old->Fill( dpel );
	hW_old->Fill( sqrt(W2recon) );

	if(WCut){
	  hhdpel_vs_vy->Fill(vy,dpel);	
	  hhdpel_vs_xtar->Fill( xtgt, dpel );
	  ep_incident_avg += T_pincident;
	  n_ep_inc++;
	}
	
	hdpel_xfp_old->Fill( xfp, dpel );
	hdpel_yfp_old->Fill( yfp, dpel );
	hdpel_xpfp_old->Fill( thfp, dpel );
	hdpel_ypfp_old->Fill( phfp, dpel );

	hdpel_xptar_old->Fill( thtgt, dpel );
	hdpel_yptar_old->Fill( phtgt, dpel );
	hdpel_ytar_old->Fill( ytgt, dpel );
	
	
	double pincident = pelastic - MeanEloss_outgoing;
	
	
	
	// cout << "Reconstructed bend angle = " << thetabend * TMath::RadToDeg() << endl;
	// cout << "p*thetabend = " << pincident * thetabend << endl;

	//Increment sums only if dpel and W are within limits:

	hpthetabend_xptar_old->Fill( thtgt, pincident*thetabend );
	
	if( dpel >= dpelmin_fit && dpel <= dpelmax_fit &&
	    W2recon >= W2min_fit && W2recon <= W2max_fit ){
	  //cout << "dpel, W = " << dpel << ", " << Wrecon << endl;

	  hp_old->Fill( precon );
	  hpthetabend_old->Fill( pincident * thetabend );

	  hpthetabend_xfp_old->Fill( xfp, pincident*thetabend );
	  hpthetabend_yfp_old->Fill( yfp, pincident*thetabend );
	  hpthetabend_xpfp_old->Fill( thfp, pincident*thetabend );
	  hpthetabend_ypfp_old->Fill( phfp, pincident*thetabend );

	  
	  hpthetabend_yptar_old->Fill( phtgt, pincident*thetabend );
	  hpthetabend_ytar_old->Fill( ytgt, pincident*thetabend );
	  
	  hxpfp_vs_xfp_old->Fill( xfp, thfp );
	  
	}
      }

      Tout->Fill();
      

  }


  ep_incident_avg /= n_ep_inc;
  
  //Let's also automate the determination of the momentum coefficients:
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);

  c1->Divide(2,2,.001,.001);
  
  TProfile *hptheta_thtgt = new TProfile( "hptheta_thtgt", "; #theta_{tgt} ; p_{incident}#theta_{bend} (GeV*rad)", 250,-0.25,0.25);
  
  
  Tout->Project("hptheta_thtgt", "ep_incident*thetabend:thtgt", "HCALcut&&BBcut");

  double lmargin = 0.16, bmargin = 0.12, rmargin = 0.12, tmargin = 0.06;
  

  hptheta_thtgt->SetMarkerStyle(20);
  c1->cd(1)->SetGrid();
  gPad->SetLeftMargin(lmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetTopMargin(tmargin);
  
  hpthetabend_xptar_old->SetTitle("");
  hpthetabend_xptar_old->SetXTitle("#theta_{tgt} (rad)");
  hpthetabend_xptar_old->SetYTitle("p_{incident}#theta_{bend} (GeV/c*rad)");
  hpthetabend_xptar_old->Draw("colz");
  hptheta_thtgt->Draw("SAME");

  
  //  hptheta_thgt->Draw("PSAME");
  
  TF1 *fitfunc = new TF1("fitfunc", "[0]*(1.0 + [1]*x)", -0.2,0.2 );

  hptheta_thtgt->Fit(fitfunc,"qS","same",thtgtmin_fit,thtgtmax_fit);

  double A_pth_fit = fitfunc->GetParameter(0);
  double B_pth_fit = fitfunc->GetParameter(1);

  TString plotexpression;
  plotexpression.Form( "%g*(1.0+%g*thtgt)/thetabend/ep_incident-1.0", A_pth_fit, B_pth_fit );
  
  Tout->Project("hdpel_new", plotexpression.Data(), "HCALcut");
  Tout->Project("hdpel_new_nocut", plotexpression.Data(), "");
  
  fitfunc = new TF1("fitfunc","pol1",-0.004,0.004);

  //hhdpel_vs_vy->Fit(fitfunc,"qR");
  
  //draw_corr(hhdpel_vs_vy,hhdpel_vs_vy);
  //cout<<fitfunc->GetParameter(0)<<endl;
  //cout<<fitfunc->GetParameter(1)<<endl;

  //W^2 = M^2 + 2M nu - Q2
  //nu = E- Eprime
  //Q2 = 2E Eprime*(1-costheta)
  //-->W^2 = M^2 + 2M * (E-Eprime) - 2E Eprime * (1-pz/p)
  plotexpression.Form("sqrt(pow(0.938272,2)+2.*.938272*(%g-(%g*(1.0+%g*thtgt)/thetabend))-2.0*%g*(%g*(1.0+%g*thtgt)/thetabend)*(1.0-cos(etheta)))",ebeam-MeanEloss, A_pth_fit, B_pth_fit, ebeam-MeanEloss, A_pth_fit, B_pth_fit);

  Tout->Project("hW_new",plotexpression.Data(), "HCALcut");
  Tout->Project("hW_new_nocut",plotexpression.Data(), "");

  c1->cd(2)->SetGrid();
  hdxdy_HCAL_old->Draw("colz");

  TEllipse E;
  E.SetFillStyle(0);
  E.SetLineColor(6);
  E.DrawEllipse(dy_p[0],dx_p[0],nsigma*dy_p[1],nsigma*dx_p[1],0,360,0);

  
  
  c1->cd(3)->SetGrid();
  gPad->SetLeftMargin(lmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetTopMargin(tmargin);

  hdpel_new_nocut->SetLineColor( kGreen+2 );
  //hdpel_new_nocut->SetFillColorAlpha( kGreen+1, 0.2 );
  
  
  hdpel_new->SetLineColor(1);
  hdpel_new->SetLineWidth(1);
  hdpel_old->SetLineColor(4);
  hdpel_old->SetLineWidth(1);
  hdpel_new->SetFillColorAlpha(1,0.7);
  hdpel_new->Draw();
  //hdpel_new->GetYaxis()->SetRangeUser(0,1.1*hdpel_new_nocut->GetMaximum());
  //hdpel_new_nocut->Draw("same");
  
  
  hdpel_new->Fit("gaus","q","",-0.04,0.03);
  hdpel_old->Draw("SAME");

  TLegend *leg = new TLegend(0.18,0.72,0.5,0.89);
  leg->SetLineColor(0);  
  leg->AddEntry(hdpel_new_nocut,"no HCAL cuts","l");
  leg->AddEntry(hdpel_old,"old coeff w/HCal cuts","l");
  leg->AddEntry(hdpel_new,"new coeff w/HCal cuts","l");
  leg->Draw("same");  

  c1->cd(4)->SetGrid();
  gPad->SetLeftMargin(lmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetTopMargin(tmargin);

  hW_new_nocut->SetLineColor( kGreen+1 );
  hW_new_nocut->SetFillColorAlpha( kGreen+1, 0.2 );
  hW_new_nocut->Draw();
  
  hW_new->SetLineColor(1);
  hW_new->SetLineWidth(1);
  hW_new->SetFillColorAlpha(1,0.7);
  hW_old->SetLineColor(4);
  hW_old->SetLineWidth(1);

  TLine L;
  L.SetLineColor(2);
  L.SetLineWidth(2);
  //hW_new->Draw("same");
  hW_new->Draw();
  L.DrawLine(0.938272,0,0.938272,hW_new->GetMaximum());
  //hW_new->GetXaxis()->SetRangeUser(0.

  // hW_new_nocut->SetLineColor(kGreen+1);
  // hW_new_nocut->Draw("SAME");
  
  hW_old->Draw("SAME");

  ///////////////// ~~~~~~~~~This section removes the y beam correlation~~~~~~~~~~~//////////////////////////

  //Now that A & B coefficients are found we will fit the beam y dependance and create coefficients for that
  TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
  c2->Divide(2,2);

  c2->cd(1);
  
  //Draw correlation and do a linear fit on the results
  hhdpel_vs_vy->Draw("colz");
  vector<double> lincuts = {0.01,-0.08,0.1,0}; //Use for GEN3
  if(cfg == "GEN2") lincuts = {0.00,-0.05,0.04,-0.01}; //Use for GEN2
  if(cfg == "GEN4b") lincuts = {-0.05,-0.05,0.04,0.04}; //Use for GEN4
  double fitresult[2];
  draw_corr(hhdpel_vs_vy,fitresult,hhdpel_vs_vy_temp,lincuts);
    
  c2->cd(2);
  
  //Since the fit above is for y = p_recon/p_elastic we must multiply the results by p_elastic so the p_recon is reconstructed
  double p_b_correct = ep_incident_avg*fitresult[0];
  double p_m_correct = ep_incident_avg*fitresult[1];
  
  //Plot the results with the new coefficients
  plotexpression.Form( "(%g*(1.0+%g*thtgt)/thetabend - (%g + %g*vy))/ep_incident-1.0:vy", A_pth_fit, B_pth_fit, p_b_correct, p_m_correct);
  Tout->Project("hhdpel_vs_vy_new", plotexpression.Data(), "HCALcut && WCut");
  hhdpel_vs_vy_new->Draw("colz");

  //For comparison look at how x_tg is affected by these new variables
  c2->cd(3);
  hhdpel_vs_xtar->Draw("colz");

  c2->cd(4);
  
  plotexpression.Form( "(%g*(1.0+%g*thtgt)/thetabend - (%g + %g*vy))/ep_incident-1.0:xtgt", A_pth_fit, B_pth_fit, p_b_correct, p_m_correct);
  Tout->Project("hhdpel_vs_xtar_new", plotexpression.Data(), "HCALcut && WCut");
  hhdpel_vs_xtar_new->Draw("colz");

  //This commented out section is if we wanted to fit x_tg instead of y beam. It does not work very well
  /*
  c2->cd(1);
    
  hhdpel_vs_xtar->Draw("colz");
  double fitresult[2];
  draw_corr(hhdpel_vs_xtar,fitresult,hhdpel_vs_xtar_temp);
  
  c2->cd(2);

  double p_b_correct = ep_incident_avg*fitresult[0];
  double p_m_correct = ep_incident_avg*fitresult[1];
  
  plotexpression.Form( "(%g*(1.0+%g*thtgt)/thetabend - (%g + %g*xtgt))/ep_incident-1.0:xtgt", A_pth_fit, B_pth_fit, p_b_correct, p_m_correct);
  
  Tout->Project("hhdpel_vs_xtar_new", plotexpression.Data(), "HCALcut && WCut");
  hhdpel_vs_xtar_new->Draw("colz");
  //draw_corr(hhdpel_vs_xtar_new,fitresult);
  
  c2->cd(3);
  hhdpel_vs_vy->Draw("colz");

  c2->cd(4);
  plotexpression.Form( "(%g*(1.0+%g*thtgt)/thetabend - (%g + %g*xtgt))/ep_incident-1.0:vy", A_pth_fit, B_pth_fit, p_b_correct, p_m_correct);
  Tout->Project("hhdpel_vs_vy_new", plotexpression.Data(), "HCALcut && WCut");
  hhdpel_vs_vy_new->Draw("colz");
  */
  /////////////////////////////////////////////////////////////////////////////////
  
  //File containing old optics coefficients: 
  ifstream foldcoeffs(fname_oldcoeffs.Data());

  TString newcoeffs_fname = outputfilename;
  newcoeffs_fname.ReplaceAll(".root",".dat");
  
  ofstream fnewcoeffs( newcoeffs_fname.Data() );

  fnewcoeffs << "bb.preconflag = 1" << endl;
  fnewcoeffs << "# NOTE on new momentum reconstruction formalism: 1st-order momentum is calculated from" << endl
	     << "#   p*thetabend = A_pth1*(1.0 + (B_pth1+C_pth1*bb.magdist)*thtgt)" << endl
	     << "# and momentum expansion coefficients are for delta = p(first order) * (1+delta)" << endl
	     << "# In most of the acceptance we don't need any corrections beyond the first-order model" << endl
	     << endl;
  

  TString currentline;
  currentline.Form("bb.A_pth1 = %15.9g", A_pth_fit );
  fnewcoeffs << currentline << endl;
  currentline.Form("bb.B_pth1 = %15.9g", B_pth_fit );
  fnewcoeffs << currentline << endl;
  fnewcoeffs << "bb.C_pth1 = 0.0" << endl << endl;
  currentline.Form("bb.A_pvy = %15.9g", p_b_correct );
  fnewcoeffs << currentline << endl;
  currentline.Form("bb.B_pvy = %15.9g", p_m_correct );
  fnewcoeffs << currentline << endl;


  fnewcoeffs << "#NOTE: angle and vertex reconstruction coefficients are from " << fname_oldcoeffs << endl << endl;
  fnewcoeffs << "bb.optics_parameters = " << endl;
  
  
  //TString currentline;
  while( currentline.ReadLine( foldcoeffs ) ){
    std::string thisline( currentline.Data() );
    std::istringstream sline( thisline );

    double coeffs[4];
    int expon[5];
    for( int i=0; i<4; i++ ){
      sline >> coeffs[i];
    }
    for( int i=0; i<5; i++ ){
      sline >> expon[i];
    }

    //Set all p*thetabend coefficients to zero: 
    coeffs[3] = 0.0; 
    TString newcoeffs;
    newcoeffs.Form(" %15.8g %15.8g %15.8g %15.8g    %d %d %d %d %d ",
		   coeffs[0], coeffs[1], coeffs[2], coeffs[3],
		   expon[0], expon[1], expon[2], expon[3], expon[4] );

    fnewcoeffs << newcoeffs << endl;
  }


  TCanvas *c5 = new TCanvas("c5","",1200,800);
  TPaveText *pt = new TPaveText(0.01,0.30,0.99,0.99,"ndc");
  //pt->AddText(runinfo);
  //if(runnums.size() >= 13) pt->AddText(runinfo2);
  //pt->AddText(cutslist);
  pt->AddText("");
  pt->AddText("");
  pt->AddText("");
  pt->AddText("New DB settings:");
  pt->AddText(Form("bb.A_pth1 = %15.9g", A_pth_fit ));
  pt->AddText(Form("bb.B_pth1 = %15.9g", B_pth_fit ));
  pt->AddText("bb.C_pth1 = 0.0");
  pt->AddText(Form("bb.A_pvy = %15.9g", p_b_correct ));
  pt->AddText(Form("bb.B_pvy = %15.9g", p_m_correct ));

  pt->SetFillColor(0);
  pt->Draw();

  TString plotsfilename = "../../plots/MomentumFit_" + cfg + ".pdf";
  
  
  c5->Print(plotsfilename + "(");
  c1->Print(plotsfilename);
  c2->Print(plotsfilename + ")");
  
  

  // elist->Delete();

  fout->Write();
}
