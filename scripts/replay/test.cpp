void test(const char *filepath){
	TFile* file = TFile::Open(filepath);
	TTree* tree = (TTree*)file->Get("T");

	double sbs_hcal_e[100];
	double bb_ps_e;
	double bb_tr_n;
	
	tree->SetBranchAddress("sbs.hcal.e",&sbs_hcal_e);
	tree->SetBranchAddress("bb.ps.e",&bb_ps_e);
	tree->SetBranchAddress("bb.tr.n",&bb_tr_n);	

	TH1D *h_sbs_hcal_e = new TH1D("h_sbs_hcal_e","sbs.hcal.e",2000,0,2);

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);
		if(sbs_hcal_e[0]>0 and bb_ps_e>0.2 and bb_tr_n>0){
		h_sbs_hcal_e->Fill(sbs_hcal_e[0]);
		}
	}
	
	TCanvas* c = new TCanvas("c","c",1600,1200);
	h_sbs_hcal_e->Draw();

}
