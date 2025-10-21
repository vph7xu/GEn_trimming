void Raw_asymmetry(const char *filename,const char *printfilename,const char *kin ){
	TFile *file = TFile::Open(filename);
	TTree *tree = (TTree*)file->Get("Tout");

 	int runnum = 0; 
	int helicity = 0;
        int IHWP = 0;


	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("IHWP",&IHWP);



}
