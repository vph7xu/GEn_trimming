
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
  double varval[maxHits];
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

  map<int,int> runs_map;

  int nruns = runnums_set.size();
  TH2D *hhvar_runnum[plots.size()][nmodules];

  for(int irun=0; irun<nruns; irun++){
    runs_map[runnums_set[irun]] = irun;
  }

  T->SetBranchStatus("*",0);
  
  // GEM hits
  double nhits, module[maxHits];
  std::vector<std::string> gemvar = {"ngoodhits","module"}; 
  std::vector<void*> gemvar_mem = {&nhits,&module};

  for(int ivar=0; ivar < plots.size(); ivar++){
    gemvar.push_back(plots[ivar].varname);
    gemvar_mem.push_back(&plots[ivar].varval);

    for( int imod=0; imod < nmodules; imod++){
      hhvar_runnum[ivar][imod] = new TH2D(Form("hh%s_runnum_m%i",plots[ivar].varname.c_str(),imod),Form(plots[ivar].title,imod),nruns,0,nruns,100,plots[ivar].yrange[0],plots[ivar].yrange[1]);
    }
  }

  setrootvar::setbranch(T, "bb.gem.hit", gemvar, gemvar_mem);

  long nevent = 0, nevents = T->GetEntries(); 
  int treenum = 0, currenttreenum = 0, currentrunnum = 0;
  int pastrunnum = runnums_set[0];  
  
  while (T->GetEntry(nevent++)) {
  
    // print progress 
    //if( nevent % 1000 == 0 ) std::cout << nevent*100.0/nevents << "% \r";
    //std::cout.flush();

    currenttreenum = T->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      string s = T->GetFile()->GetName();
      int start = s.find("_stream0");
      start -= 4;
      int end = start + 4;
      currentrunnum = stoi(s.substr(start,end - start));
    }

    for(int ihit=0; ihit < nhits; ihit++){
      for(int ivar=0; ivar < plots.size(); ivar++){
	hhvar_runnum[ivar][(int)module[ihit]]->Fill(runs_map[currentrunnum],plots[ivar].varval[ihit]*plots[ivar].scale);
      }
    }
  }
  
  TString outfile = "../../plots/" + GEN_config + "_GEM_quality_check_set" + set+ ".pdf";
  
  TCanvas *c1 = new TCanvas("c1","",1200,1200);
  TCanvas *c2 = new TCanvas("c2","",1200,1200);
  TCanvas *c3 = new TCanvas("c3","",1200,1200);
  TCanvas *c4 = new TCanvas("c4","",1200,1200);
  c1->Divide(1,2);
  c2->Divide(1,2);
  c3->Divide(1,2);
  c4->Divide(1,2);
  
  for(int ivar=0; ivar < plots.size(); ivar++){
    for(int imod=0; imod < nmodules; imod++){
      Custm2DRnumHisto(hhvar_runnum[ivar][imod],runnums_set);
      
      if(imod < 2) c1->cd(imod + 1);    
      if(imod >= 2 && imod < 4) c2->cd(imod - 1);
      if(imod >= 4 && imod < 6) c3->cd(imod - 3);    
      if(imod >= 6 && imod < 8) c4->cd(imod - 5);    
      hhvar_runnum[ivar][imod]->Draw("colz");
      
    }
   
    if(ivar == 0){
      c1->Print(outfile + "(");
      c2->Print(outfile);
      c3->Print(outfile);
      c4->Print(outfile);
    }
    else if(ivar == plots.size() - 1){
      c1->Print(outfile);
      c2->Print(outfile);
      c3->Print(outfile);
      c4->Print(outfile + ")");
    }
    else{
      c1->Print(outfile);
      c2->Print(outfile);
      c3->Print(outfile);
      c4->Print(outfile);
    }
  }

}

void GEM_quality_checks(const char *configfilename){

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

  //// U Residual //////
  plot.varname = "residu";
  plot.title = "U Residuals vs Run Number Module %i;Run Number;U Residual (#mum)";
  plot.yrange[0] = -500;
  plot.yrange[1] = 500;
  plot.scale = 1000000;
  plots.push_back(plot);

  //// V Residual //////
  plot.varname = "residv";
  plot.title = "V Residuals vs Run Number Module %i;Run Number;V Residual (#mum)";
  plot.yrange[0] = -500;
  plot.yrange[1] = 500;
  plot.scale = 1000000;
  plots.push_back(plot);

  //// ADC Asymm //////
  plot.varname = "ADCasym";
  plot.title = "ADC Asymm vs Run Number Module %i;Run Number;ADC Asymm";
  plot.yrange[0] = -1;
  plot.yrange[1] = 1;
  plot.scale = 1;
  plots.push_back(plot);

  //// ADC U Max Samples //////
  plot.varname = "ADCmaxsampU";
  plot.title = "ADC U Max Sample vs Run Number Module %i;Run Number;ADC U Max Sample";
  plot.yrange[0] = 0;
  plot.yrange[1] = 3000;
  plot.scale = 1;
  plots.push_back(plot);

  //// ADC V Max Samples //////
  plot.varname = "ADCmaxsampV";
  plot.title = "ADC V Max Sample vs Run Number Module %i;Run Number;ADC V Max Sample";
  plot.yrange[0] = 0;
  plot.yrange[1] = 3000;
  plot.scale = 1;
  plots.push_back(plot);

  //// U Time //////
  plot.varname = "UtimeMaxStrip";
  plot.title = "U Time Max Strip vs Run Number Module %i;Run Number;U Time (ns)";
  plot.yrange[0] = 0;
  plot.yrange[1] = 120;
  plot.scale = 1;
  plots.push_back(plot);

  //// V Time //////
  plot.varname = "VtimeMaxStrip";
  plot.title = "V Time Max Strip vs Run Number Module %i;Run Number;V Time (ns)";
  plot.yrange[0] = 0;
  plot.yrange[1] = 120;
  plot.scale = 1;
  plots.push_back(plot);

  //// Cluster Time delta //////
  plot.varname = "deltat";
  plot.title = "Cluster #Deltat vs Run Number Module %i;Run Number;Cluster #Deltat (ns)";
  plot.yrange[0] = -20;
  plot.yrange[1] = 20;
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
      std::string rfname = rootfile_dir + Form("/*%i*seg0_0*.root",runnums_set[i]);
      T->Add(rfname.c_str());    
    }
    if(nruns_this_set != 0) analyze_subset(T, iset + 1,runnums_set);
  }

}
 
