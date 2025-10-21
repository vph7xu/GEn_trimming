//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is to check the raster calibration.
//   The user must write down the raster limits at the beginning and
//   then it will loop over all runs and check that the raster is 
//   stable within these limits. It should be update in the future
//   to be automated.
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
#include "../../dflay/src/JSONManager.cxx"


void raster_check_calibration(const char *configfilename){

  gErrorIgnoreLevel = kError;

  //These values are taken from the SBS-replay DB
  double xmin_Raster_curr = 34720;
  double xmax_Raster_curr = 53590;
  double ymin_Raster_curr = 34705;
  double ymax_Raster_curr = 52085;

  double xmin_Raster2_curr = 43755;
  double xmax_Raster2_curr = 44111;
  double ymin_Raster2_curr = 42509;
  double ymax_Raster2_curr = 42865;

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
  
  // setting up global cuts
  std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TCut globalcut = gcut.c_str();

  TString conf = jmgr->GetValueFromKey_str("GEN_config");

  if(conf == "GEN3"){
    xmin_Raster_curr = 29106;
    xmax_Raster_curr = 59225;
    ymin_Raster_curr = 29607;
    ymax_Raster_curr = 57139;
  }

  if(conf == "GEN4"){
    xmin_Raster_curr = 34000;
    xmax_Raster_curr = 52700;
    ymin_Raster_curr = 36100;
    ymax_Raster_curr = 50700;

    xmin_Raster2_curr = 35800;
    xmax_Raster2_curr = 52000;
    ymin_Raster2_curr = 35100;
    ymax_Raster2_curr = 49900;

  }

  //Loop over all runs
  if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();
  for (int i=0; i<nruns; i++) {
    TString rfname = rootfile_dir + Form("/e1209016_fullreplay_%i_stream0_2_seg0_*.root",runnums[i]);
    //TString rfname = rootfile_dir + Form("/*%i*",runnums[i]);

    TChain *T = new TChain("T");
    T->Add(rfname);

    //Draw the raster x current
    TH1F *h = new TH1F("h","",200,200000,60000);
    T->Draw("Lrb.Raster.rawcur.x>>h","","");
    
    double xmin = h->GetBinCenter(h->FindFirstBinAbove());
    double xmax = h->GetBinCenter(h->FindLastBinAbove());

    //Draw the raster y current
    T->Draw("Lrb.Raster.rawcur.y>>h","","");
    
    double ymin = h->GetBinCenter(h->FindFirstBinAbove());
    double ymax = h->GetBinCenter(h->FindLastBinAbove());
    
    //If a run has a very weird raster tell the user
    if(abs(xmin - xmin_Raster_curr) > 200 || abs(xmax - xmax_Raster_curr) > 200 || 
       abs(ymin - ymin_Raster_curr) > 200 || abs(ymax - ymax_Raster_curr) > 200){
      cout<<"Inconsistent Raster (upstream) current found!!!"<<endl;
      cout<<"Run number "<<runnums[i]<<endl;
      cout<<"Raster x min: "<<xmin<<endl;
      cout<<"Raster x max: "<<xmax<<endl;
      cout<<"Raster y min: "<<ymin<<endl;
      cout<<"Raster y max: "<<ymax<<endl;
      cout<<endl;
    }

    //Repeat the same steps above for Raster2
    T->Draw("Lrb.Raster2.rawcur.x>>h","","");
    
    xmin = h->GetBinCenter(h->FindFirstBinAbove());
    xmax = h->GetBinCenter(h->FindLastBinAbove());

    T->Draw("Lrb.Raster2.rawcur.y>>h","","");
    
    ymin = h->GetBinCenter(h->FindFirstBinAbove());
    ymax = h->GetBinCenter(h->FindLastBinAbove());

    if(abs(xmin - xmin_Raster2_curr) > 200 || abs(xmax - xmax_Raster2_curr) > 200 || 
       abs(ymin - ymin_Raster2_curr) > 200 || abs(ymax - ymax_Raster2_curr) > 200){
      cout<<"Inconsistent Raster2 (downstream) current found!!!"<<endl;
      cout<<"Run number "<<runnums[i]<<endl;
      cout<<"Raster x min: "<<xmin<<endl;
      cout<<"Raster x max: "<<xmax<<endl;
      cout<<"Raster y min: "<<ymin<<endl;
      cout<<"Raster y max: "<<ymax<<endl;
      cout<<endl;
    }

  }



}
