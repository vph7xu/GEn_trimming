#include <iostream>
#include <TFile.h>
#include "../../include/gen-ana.h"
#include <fstream>
//#include "../../src/JSONManager.cpp"

void get_sbs_runs(string configfilename){

  string configdir = "../../config/";
  gErrorIgnoreLevel = kError; 


  JSONManager *jmgr = new JSONManager(configdir+configfilename);

  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  std::string evio_dir = "/mss/halla/sbs/GEnII/raw/";

  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  //int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana");
  int nruns = runnums.size();
  

  // Define a clock to get macro processing time
  // TStopwatch *sw = new TStopwatch(); sw->Start();
  
  // reading input config file ---------------------------------------
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(configdir + configfilename,1);
  //
  // parsing trees
  //TChain *C = LoadRawRootFiles(kin_info, 1);
  
  //const char sbsgem = "hsbs_gem_m2_nclustU_good";
  std::cout<<"number of runs : "<<nruns<<endl;
  //std::cout<<"runs : "<<runnums<<endl;
  
  string suffix = "_GEnII-3stream_runs.txt";
  string outputfilepath = configfilename+suffix;

  std::ofstream outputfile(outputfilepath);

  for (int i = 0; i<nruns; i++){  
    std::ifstream fileStream(Form("/mss/halla/sbs/GEnII/raw/e1209016_%d.evio.2.0",runnums[i]));
    
    if (fileStream.is_open()){
      std::cout<<"GEnII-3stream run"<<endl;
      outputfile<<runnums[i] <<endl;
    }
    else{
      std::cout<<"Not a GEnII-3stream run"<<endl;
    }
    fileStream.close();
    
    //..............Since mass replays excluded SBS GEMs, I cannot use the following method..................//

    //open the file
    //TFile* file = TFile::Open(Form("/cache/halla/sbs/prod/GEnII/pass1/GEN2/He3/rootfiles/e1209016_fullreplay_%d_stream0_2_seg0_0.root ",runnums[i]));
   
    //if file is not there
    //if(!file||file->IsZombie()){
      //std::cerr<<"Error : File is not there"<<std::endl;
      //continue;
    //}   
    
    //open the tree
    //TTree *T = dynamic_cast<TTree*>(file->Get("T"));
    
    //if tree is not there
    //if(!T){
      //std::cerr<<"Tree is not there in the file"<<endl;
      //file->Close();
      //continue;
    //}    
    
    //look for the branch
    //TBranch *b = T->GetBranch("sbs.tr.n");
    //if(b){
      //std::cout <<"SBS gems are there"<<endl;
      //file->Close();
      //continue;
    //}
    //else{
      //std::cout<<"SBS gems are not there"<<endl;
    //}
    //file->Close();

  }
  outputfile.close();
}
