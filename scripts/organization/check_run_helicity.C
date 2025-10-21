#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "../include/gen-ana.h"






void check_run_helicity(TString cfg_name){

  ifstream configfile("../config/" + cfg_name);
 
  TString run;
  TString Rootfiles = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEN2/He3/";

  vector<TString> run_list;
  vector<double> helicity_ratio;

  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      
      int ntokens = tokens->GetEntries();
      if( ntokens == 1 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	run = skey;

	TChain *T = new TChain("T");
	cout<<"Currently on run "<<run<<endl;

	T->Add(Rootfiles + "good_helicity/e1209016_fullreplay_*" + run + "*");
	
	run_list.push_back(run);
	helicity_ratio.push_back(T->GetEntries("scalhel.hel == 0")*1.0 / (T->GetEntries("scalhel.hel == -1") + T->GetEntries("scalhel.hel == 1")));
	
      }
      
      

    }
  }

  for( int irun=0; irun < run_list.size(); irun++)
    cout<<run_list[irun]<<" "<<helicity_ratio[irun]<<endl;
    
  cout<<"\n\n"<<"Bad Helicities:"<<endl;

  for( int irun=0; irun < run_list.size(); irun++){
    if(helicity_ratio[irun] > 0.03){ 
      cout<<run_list[irun]<<" "<<helicity_ratio[irun]<<endl;

      gSystem->Exec("mv " + Rootfiles + "good_helicity/*" + run_list[irun] + "* " + Rootfiles + "bad_helicity");
    }
  }

} 
