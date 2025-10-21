
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

int nmodules = 8;
const int maxHits=1000;

std::string GEN_config;
std::vector<int> runnums; 

struct plotinfo{
  string varname;
  TString title;
  double yrange[2];
  double varval;
  double scale;
};

vector<plotinfo> plots;


// Customizes 2D histos with run # on the X-axis
void Custm2DRnumHisto(TH2D* h, std::vector<int> const & lrnum)
{
  h->SetStats(0);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetRange(1,lrnum.size());
  
  for (int i=0; i<lrnum.size(); i++) {
    h->GetXaxis()->SetBinLabel(i+1,Form("%d",lrnum[i]));
  }
  if (lrnum.size()>15) h->LabelsOption("v", "X");
}


void analyze_subset(TChain *T, int set, std::vector<int> runnums_set){

  int nruns = runnums_set.size();
  TH2D *hhvar_runnum[plots.size()];

  T->SetBranchStatus("*",0);
  
  for(int ivar=0; ivar < plots.size(); ivar++){
    T->SetBranchStatus(plots[ivar].varname.c_str(),"1");
    T->SetBranchAddress(plots[ivar].varname.c_str(),&plots[ivar].varval);

    hhvar_runnum[ivar] = new TH2D(Form("hh%s_runnum",plots[ivar].varname.c_str()),plots[ivar].title,nruns,0,nruns,100,plots[ivar].yrange[0],plots[ivar].yrange[1]);

  }

  long nevent = 0, nevents = T->GetEntries(); 
  int treenum = 0, currenttreenum = 0, currentrunnum = 0;
  int irun = 0;
  int pastrunnum = runnums_set[0];
  
  while (T->GetEntry(nevent++)) {
  
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent*100.0/nevents << "% \r";
    std::cout.flush();

    currenttreenum = T->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      string s = T->GetFile()->GetName();
      int start = s.find("_stream0");
      start -= 4;
      int end = start + 4;
      int currentrunnum = stoi(s.substr(start,end - start));
      if(currentrunnum != pastrunnum){
	irun++;
	pastrunnum = currentrunnum;
      }
    }
    
    for(int ivar=0; ivar < plots.size(); ivar++){
      hhvar_runnum[ivar]->Fill(irun,plots[ivar].varval*plots[ivar].scale);
    }
  }
  
  TString outfile = "../../plots/" + GEN_config + "_raster_quality_check_set" + set+ ".pdf";
  
  TCanvas *c1 = new TCanvas("c1","",1200,1200);
  TCanvas *c2 = new TCanvas("c2","",1200,1200);
  c1->Divide(1,2);
  c2->Divide(1,2);
  
  for(int ivar=0; ivar < plots.size(); ivar++){
    Custm2DRnumHisto(hhvar_runnum[ivar],runnums_set);
      
    if(ivar < 2) c1->cd(ivar + 1);    
    if(ivar >= 2) c2->cd(ivar - 1);    
    
    hhvar_runnum[ivar]->Draw("colz");
      
  }
   
  c1->Print(outfile + "(");
  c2->Print(outfile + ")");

}

void raster_quality_check(const char *configfilename){

  int nruns_set = 80;  //number of runs that can fit on one pdf

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  GEN_config = jmgr->GetValueFromKey_str("GEN_config");
  jmgr->GetVectorFromKey<int>("runnums",runnums);
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
 
  // setting up global cuts
  std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TCut globalcut = gcut.c_str();
  
  if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();
  if(nruns < runnums.size()) runnums.resize(nruns);

  plotinfo plot;

  //// Raster X //////
  plot.varname = "Lrb.Raster.rawcur.x";
  plot.title = "Upstream Raster X vs Run Number;Run Number;Raster X Current";
  plot.yrange[0] = 25000;
  plot.yrange[1] = 65000;
  plot.scale = 1;
  plots.push_back(plot);

  //// Raster Y //////
  plot.varname = "Lrb.Raster.rawcur.y";
  plot.title = "Upstream Raster Y vs Run Number;Run Number;Raster Y Current";
  plot.yrange[0] = 25000;
  plot.yrange[1] = 65000;
  plot.scale = 1;
  plots.push_back(plot);

  //// Raster2 X //////
  plot.varname = "Lrb.Raster2.rawcur.x";
  plot.title = "Downstream Raster X vs Run Number;Run Number;Raster X Current";
  plot.yrange[0] = 43000;
  plot.yrange[1] = 45000;
  plot.scale = 1;
  plots.push_back(plot);

  //// Raster2 Y //////
  plot.varname = "Lrb.Raster2.rawcur.y";
  plot.title = "Downstream Raster Y vs Run Number;Run Number;Raster Y Current";
  plot.yrange[0] = 42000;
  plot.yrange[1] = 43000;
  plot.scale = 1;
  plots.push_back(plot);

  int nsets = nruns/nruns_set + 1; //Calculate number of sets
  int nremainder = nruns - nruns_set*(nsets - 1); 

  int irun = 0; 

  for(int iset=0; iset < nsets; iset++){

    int nruns_this_set = nruns_set;
    if(iset == nsets - 1) nruns_this_set = nremainder;
    
    TChain *T = new TChain("T");
    std::vector<int> runnums_set;
  
    for (int i=0; i<nruns_this_set; i++) {
      runnums_set.push_back(runnums[++irun]);
      std::string rfname = rootfile_dir + Form("/*%i*seg0_0.root",runnums_set[i]);
      T->Add(rfname.c_str()); 

      rfname = rootfile_dir + Form("/*%i*seg1_1.root",runnums_set[i]);
      T->Add(rfname.c_str()); 
    }
    if(nruns_this_set != 0) analyze_subset(T, iset + 1,runnums_set);
  }

}
 
