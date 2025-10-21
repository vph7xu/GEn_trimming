#include "../../include/gen-ana.h"
#include "../../dflay/src/JSONManager.cxx"

void check_sim_optics(){

  TString rootdir = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN2/";

  TFile *file = new TFile(rootdir + "GEN2_H2_elastic_job_0.root","read");
  TTree *T = (TTree*)file->Get("T");


  int maxNtr=1000;
  T->SetMakeClass(1);
  T->SetBranchStatus("*",0);

  //Track info
  int ntracks;
  vector<double> xfp,yfp,thfp,phfp;
  std::vector<std::string> trvar = {"ntracks","X","Y","Xp","Yp"};
  std::vector<void*> trvar_mem = {&ntracks,&xfp,&yfp,&thfp,&phfp};
  setrootvar::setbranch(T,"Earm.BBGEM.Track",trvar,trvar_mem);

  //MC info
  T->SetBranchStatus("ev",1);
  fChain->SetBranchAddress("ev", &ev_count, &b_ev);

  TString Matrix_filename = "GEN2_sim_optics.txt";
  ifstream Matrix_file(Matrix_filename);
  int row_M = 0, col_M = 9;
  Matrix_file >> row_M;
  TMatrixD M(row_M,col_M);
   
  //Read in matrix file
  for(int row=0; row<row_M; row++)
    for(int col=0; col<col_M; col++) 
      Matrix_file >> M(row,col);


  int nevent = 0;

  while(T->GetEntry(nevent++)){

    if(ntracks == 1){


      double xptgt = 0;
      double yptgt = 0;
      double ytgt = 0;
      double xtgt = 0;
      double pthetabend = 0;
      double vz = 0;
      /*
      for (int row=0; row<row_M; row++){
	xptgt += M(row,0)*pow(xfp,M(row,4))*pow(yfp,M(row,5))*pow(xpfp,M(row,6))*pow(ypfp,M(row,7))*pow(xtgt,M(row,8));
	yptgt += M(row,1)*pow(xfp,M(row,4))*pow(yfp,M(row,5))*pow(xpfp,M(row,6))*pow(ypfp,M(row,7))*pow(xtgt,M(row,8));
	ytgt += M(row,2)*pow(xfp,M(row,4))*pow(yfp,M(row,5))*pow(xpfp,M(row,6))*pow(ypfp,M(row,7))*pow(xtgt,M(row,8));
	pthetabend += M(row,3)*pow(xfp,M(row,4))*pow(yfp,M(row,5))*pow(xpfp,M(row,6))*pow(ypfp,M(row,7))*pow(xtgt,M(row,8));
	
      }
      */
      
    }
  }


}
