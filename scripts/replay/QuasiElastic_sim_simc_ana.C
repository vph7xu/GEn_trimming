//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Provakar Datta
//   Modified by Sean Jeffas, sj9ry@virginia.edu
//   Updated: finite-value protections for SIMC weights/kinematics,
//            robust normalization checks, deterministic coinCut,
//            binding of bb.tr.x/y/th/ph, and Provakar RS absolute-yield
//            normalization using the SIMC summary CSV.
//
//   The purpose of this script is to analyze simulated data for
//   SBS experiments. It is almost identical to QuasiElastic_ana.C.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>

#include "TCut.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TEllipse.h"

#include "../../include/gen-ana.h"

int QuasiElastic_sim_simc_ana(
    const std::string configfilename,
    std::string filebase="/volatile/halla/sbs/vimukthi/outfiles/SIMC/noNRad/QE_sim_simc_He3en")
{
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::fstream;
  using std::getline;
  using std::string;
  using std::stringstream;
  using std::vector;

  string configdir = "../../config/";
  gErrorIgnoreLevel = kError; // Ignores ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();

  // reading input config file ---------------------------------------
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(configdir + configfilename,0);

  // parsing trees
  TChain *C = LoadRawRootFiles(kin_info, 0);
  if (!C) {
    cerr << "ERROR: LoadRawRootFiles returned a null TChain." << endl;
    return 1;
  }

  // setting up the desired SBS configuration
  TString conf = kin_info.conf;
  TString target = kin_info.target;
  TString rootfile_dir = kin_info.rootfile_dir;
  int sbsmag = kin_info.sbsmag;
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // ----------------------------------------------------------------
  // Read SIMC normalization table
  //
  // CSV columns:
  //   1  jobid
  //   2  Nthrown
  //   3  Ntried
  //   4  genvol(MeV*sr^2)
  //   5  luminosity(ub^-1)
  //   6  ebeam(GeV)
  //   7  charge(mC)
  //   8  RndmSeed
  //   9  UsingRS
  //  10  MaxWtRS(ub/MeV/sr2)
  //  11  wtGTmaxwt
  //  12  ObsMaxWtRS
  //
  // For rejection-sampled SIMC, use Provakar's absolute-yield
  // normalization:
  //
  //                  1        MaxWtRS * genvol * luminosity
  //   weight/event = -------- --------------------------------
  //                 sum Ntried            charge
  //
  // where sum Ntried is summed over all simulation jobs included
  // in the normalization table.  The CSV charge is in mC, so it is
  // converted to C below.  The resulting event weight is therefore
  // an absolute yield per Coulomb.
  // ----------------------------------------------------------------
  struct SIMCNormInfo {
    Long64_t Nthrown = 0;
    Long64_t Ntried = 0;
    double genvol = 0.0;
    double luminosity = 0.0;
    double ebeam = 0.0;
    double charge_mC = 0.0;
    Long64_t random_seed = 0;
    int usingRS = 0;
    double maxWtRS = 0.0;
    Long64_t wtGTmaxwt = 0;
    double obsMaxWtRS = 0.0;
  };

  std::map<int, SIMCNormInfo> simc_norm;
  Long64_t total_Ntried = 0;

  TString file = rootfile_dir + "simcout/" + conf + "_" + target +
                 "_elastic_n_simc_summary.csv";

  fstream file_csv;
  file_csv.open(file.Data());

  if (!file_csv.is_open()) {
    cerr << "!!!!WARNING, ERROR: File " << file << " is not found" << endl;
    return 1;
  } else {
    cout << "DB file found: " << file << endl;
  }

  string line;
  int iline = 0;

  while (getline(file_csv, line)) {
    if (line.empty()) continue;

    iline++;
    if (iline == 1) continue; // Skip header

    stringstream ss(line);
    vector<string> cols;
    string cell;

    while (getline(ss, cell, ','))
      cols.push_back(cell);

    if (cols.size() < 12) {
      cerr << "ERROR: SIMC summary CSV line " << iline
           << " has only " << cols.size()
           << " columns; expected at least 12." << endl;
      return 1;
    }

    try {
      int jobid = std::stoi(cols[0]);

      SIMCNormInfo info;
      info.Nthrown     = std::stoll(cols[1]);
      info.Ntried      = std::stoll(cols[2]);
      info.genvol      = std::stod(cols[3]);
      info.luminosity  = std::stod(cols[4]);
      info.ebeam       = std::stod(cols[5]);
      info.charge_mC   = std::stod(cols[6]);
      info.random_seed = std::stoll(cols[7]);
      info.usingRS     = std::stoi(cols[8]);
      info.maxWtRS     = std::stod(cols[9]);
      info.wtGTmaxwt   = std::stoll(cols[10]);
      info.obsMaxWtRS  = std::stod(cols[11]);

      if (info.Ntried <= 0) {
        cerr << "ERROR: non-positive Ntried for SIMC job "
             << jobid << ": " << info.Ntried << endl;
        return 1;
      }

      if (!std::isfinite(info.genvol) ||
          !std::isfinite(info.luminosity) ||
          !std::isfinite(info.charge_mC) ||
          !std::isfinite(info.maxWtRS) ||
          !std::isfinite(info.obsMaxWtRS) ||
          info.genvol <= 0.0 ||
          info.luminosity <= 0.0 ||
          info.charge_mC <= 0.0 ||
          info.maxWtRS <= 0.0) {

        cerr << "ERROR: invalid normalization values for SIMC job "
             << jobid << endl;
        return 1;
      }

      // Provakar's weight formula uses MaxWtRS and therefore assumes
      // rejection sampling.  Stop rather than silently applying the
      // wrong normalization to a non-RS sample.
      if (info.usingRS != 1) {
        cerr << "ERROR: SIMC job " << jobid
             << " has UsingRS=" << info.usingRS
             << ".  Provakar's MaxWtRS event-weight formula is only "
             << "valid for rejection-sampled jobs." << endl;
        return 1;
      }

      if (simc_norm.count(jobid)) {
        cerr << "ERROR: duplicate jobid " << jobid
             << " in SIMC summary CSV." << endl;
        return 1;
      }

      simc_norm[jobid] = info;
      total_Ntried += info.Ntried;

      // These do not abort because older production can contain an
      // occasional overflow, but they must be visible to the analyst.
      if (info.wtGTmaxwt > 0 || info.obsMaxWtRS > info.maxWtRS) {
        cerr << "WARNING: rejection-sampling maximum was exceeded for job "
             << jobid
             << "  wtGTmaxwt=" << info.wtGTmaxwt
             << "  MaxWtRS=" << info.maxWtRS
             << "  ObsMaxWtRS=" << info.obsMaxWtRS
             << endl;
      }
    }
    catch (const std::exception &e) {
      cerr << "ERROR parsing SIMC summary CSV at line "
           << iline << "\nException: " << e.what() << endl;
      return 1;
    }
  }

  file_csv.close();

  if (simc_norm.empty() || total_Ntried <= 0) {
    cerr << "ERROR: SIMC normalization table is empty or total Ntried "
         << "is non-positive." << endl;
    return 1;
  }

  cout << "SIMC normalization jobs          : " << simc_norm.size() << endl;
  cout << "Total Ntried (sum over all jobs) : " << total_Ntried << endl;
  cout << "Using Provakar RS event weight:" << endl;
  cout << "  weight = MaxWtRS * genvol * luminosity"
       << " / (charge[C] * total_Ntried)" << endl;

  double bbtheta = sbsconf.GetBBtheta_rad();

  TVector3 BB_zaxis( sin(bbtheta), 0.0, cos(bbtheta) );
  TVector3 BB_xaxis(0,-1,0);
  TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  int model = kin_info.model;

  if (model == 0)
    cout << "Using model 0 [recon. p as indep. var.] for analysis.." << endl;
  else if (model == 1)
    cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << endl;
  else if (model == 2)
    cout << "Using model 2 [4-vector calculation] for analysis.." << endl;
  else {
    cerr << "Enter a valid model number! **!**" << endl;
    return 1;
  }

  // choosing nucleon type
  TString Ntype = kin_info.Ntype;

  // setting up global cuts
  TCut globalcut = kin_info.globalcut;
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  // setting up ROOT tree branch addresses ---------------------------
  const int maxNtr = 1000;
  C->SetBranchStatus("*",0);

  // bbcal clus var
  double eSH = 0.0, atimeSH = 0.0;
  vector<string> bbcalclvar = {"e","atimeblk"};
  vector<void*> bbcalclvar_mem = {&eSH,&atimeSH};
  setrootvar::setbranch(C,"bb.sh",bbcalclvar,bbcalclvar_mem);

  double ePS = 0.0;
  setrootvar::setbranch(C,"bb.ps","e",&ePS);

  // hcal clus var
  double eHCAL = 0.0, xHCAL = 0.0, yHCAL = 0.0;
  double rblkHCAL = 0.0, cblkHCAL = 0.0, idblkHCAL = 0.0, tdctimeHCAL = 0.0;

  vector<string> hcalclvar =
    {"e","x","y","rowblk","colblk","idblk","tdctimeblk"};

  vector<void*> hcalclvar_mem =
    {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&tdctimeHCAL};

  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  // track var
  double ntrack = 0.0;
  double p[maxNtr] = {0.0};
  double px[maxNtr] = {0.0};
  double py[maxNtr] = {0.0};
  double pz[maxNtr] = {0.0};

  double xTr[maxNtr] = {0.0};
  double yTr[maxNtr] = {0.0};
  double thTr[maxNtr] = {0.0};
  double phTr[maxNtr] = {0.0};

  double vx[maxNtr] = {0.0};
  double vy[maxNtr] = {0.0};
  double vz[maxNtr] = {0.0};

  double xtgt[maxNtr] = {0.0};
  double ytgt[maxNtr] = {0.0};
  double thtgt[maxNtr] = {0.0};
  double phtgt[maxNtr] = {0.0};

  double xfp[maxNtr] = {0.0};
  double yfp[maxNtr] = {0.0};
  double thfp[maxNtr] = {0.0};
  double phfp[maxNtr] = {0.0};

  // NOTE: x/y/th/ph were declared in the old macro but never bound.
  // Bind them here so trX/trY/trTh/trPh are not uninitialized.
  vector<string> trvar = {
    "n","p","px","py","pz",
    "x","y","th","ph",
    "vx","vy","vz",
    "tg_x","tg_y","tg_th","tg_ph",
    "r_x","r_y","r_th","r_ph"
  };

  vector<void*> trvar_mem = {
    &ntrack,&p,&px,&py,&pz,
    &xTr,&yTr,&thTr,&phTr,
    &vx,&vy,&vz,
    &xtgt,&ytgt,&thtgt,&phtgt,
    &xfp,&yfp,&thfp,&phfp
  };

  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  // MC variables
  double mc_sigma = 0.0;
  double mc_omega = 0.0;
  double mc_fnucl = 0.0;
  double mc_simc_fnucl = 0.0;
  double mc_simc_weight = 0.0;

  vector<string> mc = {"mc_sigma","mc_omega","mc_fnucl","simc_fnucl","simc_Weight"};
  vector<void*> mc_mem = {&mc_sigma,&mc_omega,&mc_fnucl,&mc_simc_fnucl,&mc_simc_weight};
  setrootvar::setbranch(C,"MC",mc,mc_mem);

  double mc_vx[maxNtr] = {0.0};
  double mc_vy[maxNtr] = {0.0};
  double mc_px[maxNtr] = {0.0};
  double mc_py[maxNtr] = {0.0};
  double mc_pz[maxNtr] = {0.0};

  vector<string> mc_vertex =
    {"sdtrack_vx","sdtrack_vy","sdtrack_momx","sdtrack_momy","sdtrack_momz"};

  vector<void*> mc_vertex_mem =
    {&mc_vx,&mc_vy,&mc_px,&mc_py,&mc_pz};

  setrootvar::setbranch(C,"MC",mc_vertex,mc_vertex_mem);

  // turning on the remaining branches used by globalcut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);
  C->SetBranchStatus("sbs.hcal.nclus", 1);

  // defining the output file
  TString outFile = Form("%s_" + sbsconf.GetSBSconf() +
                         "_sbs%dp_nucleon_%s_model%d_elastic_simc.root",
                         filebase.c_str(),
                         sbsconf.GetSBSmag(),
                         Ntype.Data(),
                         model);

  TFile *fout = new TFile(outFile.Data(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    cerr << "ERROR: could not create output file " << outFile << endl;
    return 1;
  }

  // defining histograms
  TH1F *h_W = Utilities::TH1FhW("h_W");
  TH1F *h_W_cut = Utilities::TH1FhW("h_W_cut");
  TH1F *h_W_acut = Utilities::TH1FhW("h_W_acut");

  TH1D *h_dpel =
    new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);

  TH1F *h_Q2 = Utilities::TH1FhQ2("h_Q2", conf);

  vector<double> hdx_lim = kin_info.hdx_lim;
  vector<double> hdy_lim = kin_info.hdy_lim;

  TH1F *h_dxHCAL =
    new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);",
             int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);

  TH1F *h_dyHCAL =
    new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);",
             int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);

  TH1F *h_coin_time =
    new TH1F("h_coin_time", "Coincidence time (ns)", 200, 380, 660);

  TH2F *h2_rcHCAL = Utilities::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = Utilities::TH2FdxdyHCAL("h2_dxdyHCAL");
  TH2F *h2_xyHCAL_p = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_p");
  TH2F *h2_xyHCAL_n = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_n");

  // Defining interesting ROOT tree branches
  TTree *Tout = new TTree("Tout", "");

  bool WCut = false;      Tout->Branch("WCut", &WCut, "WCut/B");
  bool pCut = false;      Tout->Branch("pCut", &pCut, "pCut/B");
  bool nCut = false;      Tout->Branch("nCut", &nCut, "nCut/B");

  double weight = 0.0;      Tout->Branch("weight", &weight, "weight/D");
  double simc_weight = 0.0; Tout->Branch("simc_weight", &simc_weight, "simc_weight/D");

  double T_mc_sigma = 0.0;  Tout->Branch("mc_sigma", &T_mc_sigma, "mc_sigma/D");
  double T_mc_omega = 0.0;  Tout->Branch("mc_omega", &T_mc_omega, "mc_omega/D");
  double T_fnucl = 0.0;     Tout->Branch("fnucl", &T_fnucl, "fnucl/D");

  bool fiduCut = false;      Tout->Branch("fiduCut", &fiduCut, "fiduCut/B");
  bool coinCut = true;       Tout->Branch("coinCut", &coinCut, "coinCut/B");

  double T_ebeam = 0.0;      Tout->Branch("ebeam", &T_ebeam, "ebeam/D");

  double T_nu = 0.0;         Tout->Branch("nu", &T_nu, "nu/D");
  double T_Q2 = 0.0;         Tout->Branch("Q2", &T_Q2, "Q2/D");
  double T_W2 = 0.0;         Tout->Branch("W2", &T_W2, "W2/D");
  double T_dpel = 0.0;       Tout->Branch("dpel", &T_dpel, "dpel/D");
  double T_ephi = 0.0;       Tout->Branch("ephi", &T_ephi, "ephi/D");
  double T_etheta = 0.0;     Tout->Branch("etheta", &T_etheta, "etheta/D");
  double T_pcentral = 0.0;   Tout->Branch("pcentral", &T_pcentral, "pcentral/D");

  double T_vz = 0.0;         Tout->Branch("vz", &T_vz, "vz/D");
  double T_vx = 0.0;         Tout->Branch("vx", &T_vx, "vx/D");
  double T_vy = 0.0;         Tout->Branch("vy", &T_vy, "vy/D");

  double T_xtgt = 0.0;       Tout->Branch("xtgt", &T_xtgt, "xtgt/D");
  double T_ytgt = 0.0;       Tout->Branch("ytgt", &T_ytgt, "ytgt/D");
  double T_thtgt = 0.0;      Tout->Branch("thtgt", &T_thtgt, "thtgt/D");
  double T_phtgt = 0.0;      Tout->Branch("phtgt", &T_phtgt, "phtgt/D");

  double T_xfp = 0.0;        Tout->Branch("xfp", &T_xfp, "xfp/D");
  double T_yfp = 0.0;        Tout->Branch("yfp", &T_yfp, "yfp/D");
  double T_thfp = 0.0;       Tout->Branch("thfp", &T_thfp, "thfp/D");
  double T_phfp = 0.0;       Tout->Branch("phfp", &T_phfp, "phfp/D");

  double T_trP = 0.0;        Tout->Branch("trP", &T_trP, "trP/D");
  double T_trX = 0.0;        Tout->Branch("trX", &T_trX, "trX/D");
  double T_trY = 0.0;        Tout->Branch("trY", &T_trY, "trY/D");
  double T_trTh = 0.0;       Tout->Branch("trTh", &T_trTh, "trTh/D");
  double T_trPh = 0.0;       Tout->Branch("trPh", &T_trPh, "trPh/D");

  double T_ePS = 0.0;        Tout->Branch("ePS", &T_ePS, "ePS/D");
  double T_eSH = 0.0;        Tout->Branch("eSH", &T_eSH, "eSH/D");

  double T_eHCAL = 0.0;      Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D");
  double T_xHCAL = 0.0;      Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D");
  double T_yHCAL = 0.0;      Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D");

  double T_xHCAL_exp = 0.0;  Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D");
  double T_yHCAL_exp = 0.0;  Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D");

  double T_dx = 0.0;         Tout->Branch("dx", &T_dx, "dx/D");
  double T_dy = 0.0;         Tout->Branch("dy", &T_dy, "dy/D");
  double T_theta_pq = 0.0;   Tout->Branch("theta_pq", &T_theta_pq, "theta_pq/D");

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

  vector<double> hcal_active_area =
    cut::hcal_active_area_data();

  vector<double> hcal_safety_margin =
    cut::hcal_safety_margin(dx_p[1], dx_n[1], dy_p[1], hcal_active_area);

  // elastic cut limits
  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;

  // construct axes of HCAL CoS in Hall CoS
  double hcal_voffset = kin_info.hcal_voffset;
  double hcal_hoffset = kin_info.hcal_hoffset;

  vector<TVector3> HCAL_axes;
  kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);

  TVector3 HCAL_origin =
    sbsconf.GetHCALdist()*HCAL_axes[2] +
    hcal_voffset*HCAL_axes[0] +
    hcal_hoffset*HCAL_axes[1];

  // ----------------------------------------------------------------
  // Diagnostic counters
  // ----------------------------------------------------------------
  Long64_t n_bad_jobid      = 0;
  Long64_t n_bad_norm       = 0;
  Long64_t n_bad_weight     = 0;
  Long64_t n_bad_track      = 0;
  Long64_t n_bad_kinematics = 0;
  Long64_t n_bad_hcal       = 0;
  Long64_t n_written        = 0;

  // looping through the tree ---------------------------------------
  cout << endl;

  Long64_t nevent = 0;
  Long64_t nevents = C->GetEntries();

  int treenum = -1;
  int currenttreenum = -1;
  int run_inc = -1;

  cout << "Processing " << nevents << " events" << endl;

  while (C->GetEntry(nevent++)) {

    if (nevent % 1000 == 0)
      cout << nevent*100.0/nevents << "% \r";
    cout.flush();

    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();

    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();

      string s = C->GetFile()->GetName();

      std::size_t start = s.find("job_");
      std::size_t end = s.rfind(".root");

      if (start == string::npos || end == string::npos) {
        cerr << "\nWARNING: could not extract job number from file: "
             << s << endl;
        run_inc = -1;
      } else {
        start += 4;
        try {
          run_inc = std::stoi(s.substr(start, end-start));
        }
        catch (const std::exception &e) {
          cerr << "\nWARNING: failed to parse job number from file: "
               << s << "\n" << e.what() << endl;
          run_inc = -1;
        }
      }
    }

    bool passedgCut = GlobalCut->EvalInstance(0) != 0;
    if (!passedgCut) continue;

    // --------------------------------------------------------------
    // Look up this simulation job in the normalization table
    // --------------------------------------------------------------
    auto norm_it = simc_norm.find(run_inc);

    if (run_inc < 0 || norm_it == simc_norm.end()) {
      n_bad_jobid++;

      if (n_bad_jobid <= 10) {
        cerr << "\nWARNING: no SIMC normalization entry for run_inc/jobid = "
             << run_inc << endl;
      }

      continue;
    }

    const SIMCNormInfo &norm = norm_it->second;

    // --------------------------------------------------------------
    // MC event weighting
    // --------------------------------------------------------------
    // For rejection-sampled SIMC, the accepted events already carry
    // the generated physics distribution through the accept/reject step.
    // Do NOT multiply accepted events by MC.simc_Weight again.
    //
    // Provakar's absolute-yield normalization (slide "Absolute Yield
    // from MC") assigns the same normalization weight to every accepted
    // event from a given simulation job:
    //
    //                  MaxWtRS_job * genvol_job * luminosity_job
    //   weight/event = -----------------------------------------------
    //                    total_Ntried * charge_job[C]
    //
    // total_Ntried is the SUM of Ntried over all jobs in the CSV.

    T_mc_sigma = std::isfinite(mc_sigma) ? mc_sigma : 0.0;
    T_mc_omega = std::isfinite(mc_omega) ? mc_omega : 0.0;
    T_fnucl    = std::isfinite(mc_simc_fnucl) ? mc_simc_fnucl : 0.0;

    const double charge_C = norm.charge_mC * 1.0e-3;
    const double total_Ntried_d = static_cast<double>(total_Ntried);

    if (!std::isfinite(norm.maxWtRS) ||
        !std::isfinite(norm.genvol) ||
        !std::isfinite(norm.luminosity) ||
        !std::isfinite(charge_C) ||
        !std::isfinite(total_Ntried_d) ||
        norm.maxWtRS <= 0.0 ||
        norm.genvol <= 0.0 ||
        norm.luminosity <= 0.0 ||
        charge_C <= 0.0 ||
        total_Ntried <= 0) {

      n_bad_norm++;

      if (n_bad_norm <= 10) {
        cerr << "\nWARNING: invalid Provakar SIMC normalization"
             << "\n  jobid          = " << run_inc
             << "\n  MaxWtRS        = " << norm.maxWtRS
             << "\n  genvol         = " << norm.genvol
             << "\n  luminosity     = " << norm.luminosity
             << "\n  charge(mC)     = " << norm.charge_mC
             << "\n  charge(C)      = " << charge_C
             << "\n  total_Ntried   = " << total_Ntried
             << endl;
      }

      continue;
    }

    simc_weight =
      norm.maxWtRS *
      norm.genvol *
      norm.luminosity /
      (charge_C * total_Ntried_d);

    if (!std::isfinite(simc_weight)) {
      n_bad_weight++;

      if (n_bad_weight <= 10) {
        cerr << "\nWARNING: non-finite Provakar event weight"
             << "\n  jobid          = " << run_inc
             << "\n  MaxWtRS        = " << norm.maxWtRS
             << "\n  genvol         = " << norm.genvol
             << "\n  luminosity     = " << norm.luminosity
             << "\n  charge(C)      = " << charge_C
             << "\n  total_Ntried   = " << total_Ntried
             << "\n  event weight   = " << simc_weight
             << endl;
      }

      continue;
    }

    // Store the absolute-yield-per-Coulomb event normalization in both
    // branches used by the downstream analysis.
    weight = simc_weight;

    // --------------------------------------------------------------
    // Validate reconstructed track BEFORE using p[0], etc.
    // --------------------------------------------------------------
    if (!(ntrack > 0) ||
        !std::isfinite(p[0]) ||
        !std::isfinite(px[0]) ||
        !std::isfinite(py[0]) ||
        !std::isfinite(pz[0]) ||
        !std::isfinite(vz[0]) ||
        p[0] <= 0.0) {

      n_bad_track++;
      continue;
    }

    // kinematic parameters
    const double ebeam = sbsconf.GetEbeam();
    const double ebeam_corr = ebeam;
    const double precon = p[0];

    if (!std::isfinite(ebeam_corr) ||
        !std::isfinite(precon) ||
        precon <= 0.0) {

      n_bad_kinematics++;
      continue;
    }

    // Since precon == p[0] here, pscale is exactly 1.
    // Keep this form so an outgoing-energy-loss correction can be
    // restored later without changing the four-vector construction.
    const double pscale = precon / p[0];

    if (!std::isfinite(pscale)) {
      n_bad_kinematics++;
      continue;
    }

    // constructing the 4 vectors
    TVector3 vertex(0, 0, vz[0]);

    TLorentzVector Pe(
      0.0,
      0.0,
      ebeam_corr,
      ebeam_corr
    );

    TLorentzVector Peprime(
      px[0] * pscale,
      py[0] * pscale,
      pz[0] * pscale,
      precon
    );

    TLorentzVector PN;
    kine::SetPN(Ntype.Data(), PN);

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pcentral = kine::pcentral(ebeam_corr, etheta, Ntype.Data());

    if (!std::isfinite(etheta) ||
        !std::isfinite(ephi) ||
        !std::isfinite(pcentral) ||
        pcentral <= 0.0) {

      n_bad_kinematics++;
      continue;
    }

    double nu = 0.0;
    double pN_expect = 0.0;
    double thetaN_expect = 0.0;
    double phiN_expect = ephi + constant::pi;

    TVector3 pNhat;
    double Q2recon = 0.0;
    double W2recon = 0.0;

    if (model == 0) {

      nu = Pe.E() - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype.Data());

      if (!std::isfinite(pN_expect) || pN_expect <= 0.0) {
        n_bad_kinematics++;
        continue;
      }

      double acosterm = (Pe.E() - Peprime.Pz()) / pN_expect;
      acosterm = std::max(-1.0, std::min(1.0, acosterm));

      thetaN_expect = acos(acosterm);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);

      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());

    } else if (model == 1) {

      nu = Pe.E() - pcentral;
      pN_expect = kine::pN_expect(nu, Ntype.Data());

      if (!std::isfinite(pN_expect) || pN_expect <= 0.0) {
        n_bad_kinematics++;
        continue;
      }

      double acosterm =
        (Pe.E() - pcentral*cos(etheta)) / pN_expect;

      acosterm = std::max(-1.0, std::min(1.0, acosterm));

      thetaN_expect = acos(acosterm);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);

      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());

    } else if (model == 2) {

      TLorentzVector q = Pe - Peprime;
      nu = q.E();

      if (!std::isfinite(nu) || q.Vect().Mag2() <= 0.0) {
        n_bad_kinematics++;
        continue;
      }

      pNhat = q.Vect().Unit();
      Q2recon = -q.M2();
      W2recon = (PN + q).M2();
    }

    if (!std::isfinite(nu) ||
        !std::isfinite(Q2recon) ||
        !std::isfinite(W2recon) ||
        !std::isfinite(pNhat.X()) ||
        !std::isfinite(pNhat.Y()) ||
        !std::isfinite(pNhat.Z())) {

      n_bad_kinematics++;
      continue;
    }

    h_Q2->Fill(Q2recon, weight);

    double Wrecon = sqrt(std::max(0.0, W2recon));
    double dpel = Peprime.E()/pcentral - 1.0;

    if (!std::isfinite(Wrecon) || !std::isfinite(dpel)) {
      n_bad_kinematics++;
      continue;
    }

    h_dpel->Fill(dpel, weight);

    // MC true momentum diagnostic vectors
    if (!std::isfinite(mc_px[0]) ||
        !std::isfinite(mc_py[0]) ||
        !std::isfinite(mc_pz[0])) {
      n_bad_kinematics++;
      continue;
    }

    TVector3 pvect_true(mc_px[0],mc_py[0],mc_pz[0]);
    TVector3 pvect_BB(
      pvect_true.Dot(BB_xaxis),
      pvect_true.Dot(BB_yaxis),
      pvect_true.Dot(BB_zaxis)
    );
    (void)pvect_BB; // retained for future diagnostics

    // --------------------------------------------------------------
    // Copy reconstructed variables to output
    // --------------------------------------------------------------
    T_ebeam = Pe.E();

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pcentral = pcentral;

    T_vz = vz[0];
    T_vx = std::isfinite(mc_vx[0]) ? mc_vx[0] : 0.0;
    T_vy = std::isfinite(mc_vy[0]) ? mc_vy[0] : 0.0;

    T_xtgt = xtgt[0];
    T_ytgt = ytgt[0];
    T_thtgt = thtgt[0];
    T_phtgt = phtgt[0];

    T_xfp = xfp[0];
    T_yfp = yfp[0];
    T_thfp = thfp[0];
    T_phfp = phfp[0];

    T_trP = p[0];
    T_trX = xTr[0];
    T_trY = yTr[0];
    T_trTh = thTr[0];
    T_trPh = phTr[0];

    T_ePS = ePS;
    T_eSH = eSH;

    T_eHCAL = eHCAL;
    T_xHCAL = xHCAL;
    T_yHCAL = yHCAL;

    // Reject non-finite detector/reconstructed quantities before
    // calculating expected HCAL position or writing the event.
    if (!std::isfinite(T_xtgt) ||
        !std::isfinite(T_ytgt) ||
        !std::isfinite(T_thtgt) ||
        !std::isfinite(T_phtgt) ||
        !std::isfinite(T_xfp) ||
        !std::isfinite(T_yfp) ||
        !std::isfinite(T_thfp) ||
        !std::isfinite(T_phfp) ||
        !std::isfinite(T_trP) ||
        !std::isfinite(T_trX) ||
        !std::isfinite(T_trY) ||
        !std::isfinite(T_trTh) ||
        !std::isfinite(T_trPh) ||
        !std::isfinite(T_ePS) ||
        !std::isfinite(T_eSH) ||
        !std::isfinite(T_eHCAL) ||
        !std::isfinite(T_xHCAL) ||
        !std::isfinite(T_yHCAL)) {

      n_bad_hcal++;
      continue;
    }

    double theta_pq =
      kine::theta_pq(vertex, HCAL_origin, pNhat, xHCAL, yHCAL);

    vector<double> xyHCAL_exp;
    kine::GetxyHCALexpect(
      vertex,
      pNhat,
      HCAL_origin,
      HCAL_axes,
      xyHCAL_exp
    );

    if (xyHCAL_exp.size() < 2 ||
        !std::isfinite(theta_pq) ||
        !std::isfinite(xyHCAL_exp[0]) ||
        !std::isfinite(xyHCAL_exp[1])) {

      n_bad_hcal++;
      continue;
    }

    double dx = xHCAL - xyHCAL_exp[0];
    double dy = yHCAL - xyHCAL_exp[1];

    if (!std::isfinite(dx) || !std::isfinite(dy)) {
      n_bad_hcal++;
      continue;
    }

    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_dx = dx;
    T_dy = dy;
    T_theta_pq = theta_pq;

    // HCAL active area and safety margin cuts
    bool AR_cut =
      cut::inHCAL_activeA(xHCAL, yHCAL, hcal_active_area);

    bool FR_cut =
      cut::inHCAL_fiducial(
        xyHCAL_exp[0],
        xyHCAL_exp[1],
        sbs_kick,
        hcal_safety_margin
      );

    // HCAL cuts
    pCut =
      pow((dx-dx_p[0]) / (dx_p[1]*Nsigma_cut_dx_p), 2) +
      pow((dy-dy_p[0]) / (dy_p[1]*Nsigma_cut_dy_p), 2)
      <= 1.0;

    nCut =
      pow((dx-dx_n[0]) / (dx_n[1]*Nsigma_cut_dx_n), 2) +
      pow((dy-dy_n[0]) / (dy_n[1]*Nsigma_cut_dy_n), 2)
      <= 1.0;

    fiduCut = AR_cut && FR_cut;
    WCut = W2recon > W2min && W2recon < W2max;

    // No coincidence-time quantity is currently read for this SIMC macro.
    // Make the behavior explicit and deterministic instead of leaving
    // coinCut uninitialized.
    coinCut = true;

    // --------------------------------------------------------------
    // Final guard: do not write any event containing NaN/Inf
    // --------------------------------------------------------------
    const bool output_is_finite =
      std::isfinite(weight) &&
      std::isfinite(simc_weight) &&
      std::isfinite(T_mc_sigma) &&
      std::isfinite(T_mc_omega) &&
      std::isfinite(T_fnucl) &&
      std::isfinite(T_ebeam) &&
      std::isfinite(T_nu) &&
      std::isfinite(T_Q2) &&
      std::isfinite(T_W2) &&
      std::isfinite(T_dpel) &&
      std::isfinite(T_ephi) &&
      std::isfinite(T_etheta) &&
      std::isfinite(T_pcentral) &&
      std::isfinite(T_vz) &&
      std::isfinite(T_vx) &&
      std::isfinite(T_vy) &&
      std::isfinite(T_xtgt) &&
      std::isfinite(T_ytgt) &&
      std::isfinite(T_thtgt) &&
      std::isfinite(T_phtgt) &&
      std::isfinite(T_xfp) &&
      std::isfinite(T_yfp) &&
      std::isfinite(T_thfp) &&
      std::isfinite(T_phfp) &&
      std::isfinite(T_trP) &&
      std::isfinite(T_trX) &&
      std::isfinite(T_trY) &&
      std::isfinite(T_trTh) &&
      std::isfinite(T_trPh) &&
      std::isfinite(T_ePS) &&
      std::isfinite(T_eSH) &&
      std::isfinite(T_eHCAL) &&
      std::isfinite(T_xHCAL) &&
      std::isfinite(T_yHCAL) &&
      std::isfinite(T_xHCAL_exp) &&
      std::isfinite(T_yHCAL_exp) &&
      std::isfinite(T_dx) &&
      std::isfinite(T_dy) &&
      std::isfinite(T_theta_pq);

    if (!output_is_finite) {
      n_bad_kinematics++;
      continue;
    }

    // W cut and coincidence cut
    if (WCut && coinCut) {
      h_dxHCAL->Fill(dx,weight);
      h_dyHCAL->Fill(dy,weight);

      h2_rcHCAL->Fill(cblkHCAL, rblkHCAL, weight);
      h2_dxdyHCAL->Fill(dy, dx, weight);

      if (pCut)
        h2_xyHCAL_p->Fill(
          xyHCAL_exp[1],
          xyHCAL_exp[0] - sbs_kick,
          weight
        );

      if (nCut)
        h2_xyHCAL_n->Fill(
          xyHCAL_exp[1],
          xyHCAL_exp[0],
          weight
        );
    }

    h_W->Fill(Wrecon, weight);

    if (pCut || nCut)
      h_W_cut->Fill(Wrecon, weight);
    else
      h_W_acut->Fill(Wrecon, weight);

    Tout->Fill();
    n_written++;

  } // event loop

  cout << endl << endl;

  // ----------------------------------------------------------------
  // Summary of skipped pathological events
  // ----------------------------------------------------------------
  cout << "========== SIMC skim diagnostics ==========" << endl;
  cout << "Input entries                  : " << nevents << endl;
  cout << "Output entries written         : " << n_written << endl;
  cout << "Invalid job-index events       : " << n_bad_jobid << endl;
  cout << "Invalid normalization events   : " << n_bad_norm << endl;
  cout << "Non-finite weight events       : " << n_bad_weight << endl;
  cout << "Invalid track events           : " << n_bad_track << endl;
  cout << "Invalid kinematic events       : " << n_bad_kinematics << endl;
  cout << "Invalid HCAL/output events     : " << n_bad_hcal << endl;
  cout << "===========================================" << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
  c1->Divide(2,2);

  c1->cd(1);
  h2_dxdyHCAL->Draw("colz");

  TEllipse Ep_p;
  Ep_p.SetFillStyle(0);
  Ep_p.SetLineColor(2);
  Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(
    dy_p[0],
    dx_p[0],
    Nsigma_cut_dy_p*dy_p[1],
    Nsigma_cut_dx_p*dx_p[1],
    0,360,0
  );

  TEllipse Ep_n;
  Ep_n.SetFillStyle(0);
  Ep_n.SetLineColor(3);
  Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(
    dy_n[0],
    dx_n[0],
    Nsigma_cut_dy_n*dy_n[1],
    Nsigma_cut_dx_n*dx_n[1],
    0,360,0
  );

  c1->cd(2);
  h_W->SetLineColor(1);
  h_W->Draw();

  h_W_cut->SetLineColor(2);
  h_W_cut->Draw("same");

  h_W_acut->Draw("same");

  c1->cd(3);
  h2_xyHCAL_p->Draw("colz");
  Utilities::DrawArea(hcal_active_area);
  Utilities::DrawArea(hcal_safety_margin,4);

  c1->cd(4);
  h2_xyHCAL_n->Draw("colz");
  Utilities::DrawArea(hcal_active_area);
  Utilities::DrawArea(hcal_safety_margin,4);

  cout << "------" << endl;
  cout << " Output file : " << outFile << endl;
  cout << "------" << endl << endl;

  sw->Stop();

  cout << "CPU time elapsed = " << sw->CpuTime()
       << " s. Real time = " << sw->RealTime()
       << " s." << endl << endl;

  fout->cd();
  c1->Write();
  fout->Write();
  fout->Close();

  delete GlobalCut;
  sw->Delete();

  return 0;
}

