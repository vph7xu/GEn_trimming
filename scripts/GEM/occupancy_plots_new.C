
const int nQ2 = 3;
const int nlayers = 5;
const int nmodules = 8;

double Q2_list[nQ2] = {3.0,6.8,9.8};
TString kin_list[nQ2] = {"GEN2","GEN3","GEN4"};
int runs[nQ2] = {2164,2566,3854};
double curr_array[nQ2] = {45,45,45};
double GEM_eff[nQ2][nlayers];

int NStripsU[nlayers] = {3840, 3840, 3840, 3840, 5120};
int NStripsV[nlayers] = {3840, 3840, 3840, 3840, 1536*4};


double modx[nQ2][nmodules] = {
  {0,    -0.00187117,    -0.00744441,    -0.00902334,      -0.764072,      -0.252332,       0.258671,        0.77031},
  {0,    -0.00193194,    -0.00754019,    -0.00913784,      -0.764702,      -0.252366,         0.2587,       0.770636},
  {0,    -0.00196804,    -0.00766061,    -0.00922111,      -0.764362,       -0.25226,       0.258809,       0.770737}
};

double mody[nQ2][nmodules] = {
  {0,     -0.0112376,     0.00605578,     0.00519494,     0.00181073,    -0.00241956,    -0.00826899,     -0.0118803},
  {0,     -0.0114293,     0.00587946,     0.00491722,      0.0020337,    -0.00228423,    -0.00813535,      -0.011822},
  {0     -0.0114109,     0.00524223,     0.00493975,     0.00196191,     -0.0023144,    -0.00833631,     -0.0117745}
};

double modz[nQ2][nmodules] = {
  {0,       0.164888,       0.278845,       0.457217,        1.57534,        1.53481,        1.57823,        1.53487},
  {0,       0.165205,       0.279033,       0.457519,        1.57894,         1.5359,        1.57844,        1.53718},
  {0,       0.165557,       0.279402,       0.458364,        1.58039,        1.53848,         1.5836,        1.53889}
};

double modax[nQ2][nmodules] = {
  {0,   -0.000480217,    0.000238544,    0.000101212,    0.000382454,   -0.000134208,     0.00359775,   -0.000732322},
  {0,     0.00135372,     0.00228606,     0.00377005,    -0.00366946,    -0.00400137,    -0.00324828,    -0.00356756},
  {0,     0.00280863,     0.00316803,     0.00376198,   -0.000215649,    -0.00411266,     0.00529021,    -0.00361049}
};

double moday[nQ2][nmodules] = {
  {0,     0.00104002,     0.00142587,     0.00191641,    -0.00514573,   -0.000638233,     0.00238349,    0.000605008},
  {0,    0.000872406,     0.00128448,     0.00176611,    -0.00517063,   -0.000726478,    -0.00319687,    -0.00121113},
  {0,     0.00107349,     0.00202562,     0.00237895,    -0.00563243,    -0.00288393,    -0.00665864,   -0.000273707}
};

double modaz[nQ2][nmodules] = {
  {0,    0.000989535,    -0.00527547,    -0.00820699,    -0.00539627,    -0.00499825,    -0.00453188,    -0.00614182},
  {0,     0.00107533,     -0.0051902,     -0.0079445,    -0.00604142,    -0.00507992,    -0.00437341,    -0.00485518},
  {0,     0.00108987,    -0.00524883,    -0.00793623,     -0.0055565,    -0.00520223,    -0.00458386,     -0.0048047}
};

double mod_sizex[nmodules] = {1.5, 1.5, 1.5, 1.5, 0.512, 0.512, 0.512, 0.512};
double mod_sizey[nmodules] = {0.4, 0.4, 0.4, 0.4, 0.6144, 0.6144, 0.6144, 0.6144};


TVector3 fXax[nQ2], fYax[nQ2], fZax[nQ2];


double probability(int success, int iQ2){
  
  double result = 0;
  vector<vector<int>> hit;
  bool first_time = true;

  if(success == 5){
    result = GEM_eff[iQ2][0] * GEM_eff[iQ2][1] * GEM_eff[iQ2][2] * GEM_eff[iQ2][3] * GEM_eff[iQ2][4];
  }
  else {
    for (int i = 0; i < nlayers; ++i) {
      for (int j = 0; j < nlayers; ++j) {
	if (j != i) {
	  for (int k = 0; k < nlayers; ++k) {
	    if (k != i && k != j) {

	      if(success == 3){
		first_time = true;
		vector<int> temp_arr = {i,j,k};
		sort(temp_arr.begin(),temp_arr.end());
	      
		for(int ihit=0; ihit < hit.size();ihit++){
		  if(hit[ihit] == temp_arr) first_time = false;
		}
		hit.push_back(temp_arr);
	      }

	      for (int l = 0; l < nlayers; ++l) {
		if (l != i && l != j && l != k) {

		  if(success == 4){
		    first_time = true;
		    vector<int> temp_arr = {i,j,k,l};
		    sort(temp_arr.begin(),temp_arr.end());
	      
		    for(int ihit=0; ihit < hit.size();ihit++){
		      if(hit[ihit] == temp_arr) first_time = false;
		    }
		    hit.push_back(temp_arr);
		  }

		  for (int m = 0; m < nlayers; ++m) {
		    if (m != i && m != j && m != k && m != l) {
		    
		      vector<int> temp_arr2;
		      if(success ==3){
			temp_arr2.push_back(i); temp_arr2.push_back(j); temp_arr2.push_back(k);
		      }
		      if(success ==4){
			temp_arr2.push_back(i); temp_arr2.push_back(j); temp_arr2.push_back(k); temp_arr2.push_back(l);
		      }
		      sort(temp_arr2.begin(),temp_arr2.end());
		    
		      bool skip_this = false;
		    
		      if(!first_time){
			for(int ihit=0; ihit < hit.size();ihit++){
			  if(hit[ihit] == temp_arr2 && !first_time){ 
			    skip_this = true;
			  }
			}
		      }
		      else first_time = false;

		      if(!skip_this){
			if (success == 3) result += GEM_eff[iQ2][i] * GEM_eff[iQ2][j] * GEM_eff[iQ2][k] * (1 - GEM_eff[iQ2][l]) * (1 - GEM_eff[iQ2][m]);
			if (success == 4) result += GEM_eff[iQ2][i] * GEM_eff[iQ2][j] * GEM_eff[iQ2][k] * GEM_eff[iQ2][l] * (1 - GEM_eff[iQ2][m]);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  
  return result;
}

TVector3 TrackIntersect( int module, TVector3 track_origin, TVector3 track_direction, int iQ2){
  TVector3 modpos(modx[iQ2][module],mody[iQ2][module],modz[iQ2][module]);
  TVector3 modzaxis = fZax[iQ2];
  
  double sintersect = modzaxis.Dot( modpos - track_origin )/modzaxis.Dot( track_direction );
  
  return track_origin + sintersect * track_direction;
}

bool IsInActiveArea(int module, TVector3 point, int iQ2){

  TVector3 fOrigin(modx[iQ2][module], mody[iQ2][module], modz[iQ2][module]);
  TVector3 v = point - fOrigin;
  TVector3 det(v.Dot(fXax[iQ2]), v.Dot(fYax[iQ2]), v.Dot(fZax[iQ2]));
  

  return (abs(det.X()) <= 0.5*mod_sizex[module] && abs(det.Y()) <= 0.5*mod_sizey[module]);

}


void occupancy_plots_new(){

  double U_occu[nQ2][nlayers];
  double U_err[nQ2][nlayers];
  double efficiency[nQ2][nlayers];


  TString rootdir = "/cache/halla/sbs/prod/GEnII/pass1/";

  for(int iQ2=0; iQ2 < nQ2; iQ2++){
    
    //Initialize modules coordinates (used for efficiency calculation)
    for(int imod=0; imod < nmodules; imod++){
    
      TRotation RotTemp;

      RotTemp.RotateX( modax[iQ2][imod] * TMath::DegToRad() );
      RotTemp.RotateY( moday[iQ2][imod] * TMath::DegToRad() );
      RotTemp.RotateZ( modaz[iQ2][imod] * TMath::DegToRad() );

      fXax[iQ2].SetXYZ( RotTemp.XX(), RotTemp.YX(), RotTemp.ZX() );
      fYax[iQ2].SetXYZ( RotTemp.XY(), RotTemp.YY(), RotTemp.ZY() );
      fZax[iQ2].SetXYZ( RotTemp.XZ(), RotTemp.YZ(), RotTemp.ZZ() );
    
    }

    TH1F *hist[nlayers];
    for(int ilayer=0; ilayer < nlayers; ilayer++){
      hist[ilayer] = new TH1F(Form("hist%i",ilayer),"",100,0,1000);
    }

    TChain *C = new TChain("T");
    TChain *Tscaler = new TChain("TSsbs");

    for(int iseg=0;iseg < 100;iseg++){
      C->Add(Form(rootdir + kin_list[iQ2] + "/He3/rootfiles/*%i*seg%i_%i*",runs[iQ2],iseg,iseg));    
      Tscaler->Add(Form(rootdir + kin_list[iQ2] + "/He3/rootfiles/*%i*seg%i_%i*",runs[iQ2],iseg,iseg));    
    }

    C->SetBranchStatus("*",0);
    Tscaler->SetBranchStatus("*",0);
    
    int maxtr = 1000;
    int maxhit = 10000;
    double tr_n, hit_n;
    double tr_x[maxtr], tr_y[maxtr], tr_xp[maxtr], tr_yp[maxtr], tr_vz[maxtr], px[maxtr], py[maxtr], pz[maxtr], p[maxtr];
    double hit_tr_i[maxhit], hit_mod[maxhit], nstrips[maxhit];
    double ps_e, beam_evnum, evnum, beam_curr;

    C->SetBranchStatus("bb.tr.n",1);
    C->SetBranchStatus("bb.tr.x",1);
    C->SetBranchStatus("bb.tr.y",1);
    C->SetBranchStatus("bb.tr.th",1);
    C->SetBranchStatus("bb.tr.ph",1);
    C->SetBranchStatus("bb.tr.px",1);
    C->SetBranchStatus("bb.tr.py",1);
    C->SetBranchStatus("bb.tr.pz",1);
    C->SetBranchStatus("bb.tr.p",1);
    C->SetBranchStatus("bb.tr.vz",1);
    C->SetBranchStatus("bb.gem.nstripsu_layer",1);
    C->SetBranchStatus("bb.gem.hit.ngoodhits",1);
    C->SetBranchStatus("bb.gem.hit.trackindex",1);
    C->SetBranchStatus("bb.gem.hit.module",1);
    C->SetBranchStatus("bb.ps.e",1);
    C->SetBranchStatus("g.evnum",1);
    Tscaler->SetBranchStatus("evNumber",1);
    Tscaler->SetBranchStatus("sbs.bcm.u1.current",1);

    C->SetBranchAddress("bb.tr.n",&tr_n);
    C->SetBranchAddress("bb.tr.x",tr_x);
    C->SetBranchAddress("bb.tr.y",tr_y);
    C->SetBranchAddress("bb.tr.th",tr_xp);
    C->SetBranchAddress("bb.tr.ph",tr_yp);
    C->SetBranchAddress("bb.tr.px",px);
    C->SetBranchAddress("bb.tr.py",py);
    C->SetBranchAddress("bb.tr.pz",pz);
    C->SetBranchAddress("bb.tr.p",p);
    C->SetBranchAddress("bb.tr.vz",tr_vz);
    C->SetBranchAddress("bb.gem.nstripsu_layer",nstrips);
    C->SetBranchAddress("bb.gem.hit.ngoodhits",&hit_n);
    C->SetBranchAddress("bb.gem.hit.trackindex",hit_tr_i);
    C->SetBranchAddress("bb.gem.hit.module",hit_mod);
    C->SetBranchAddress("bb.ps.e",&ps_e);
    C->SetBranchAddress("g.evnum",&evnum);
    Tscaler->SetBranchAddress("evNumber",&beam_evnum);
    Tscaler->SetBranchAddress("sbs.bcm.u1.current",&beam_curr);

    int set_max_events = 10000000;
    int nevent_scaler = 0;
    int nevent = 0;
    int nevent_start = 0;
    int nevent_end = -1;
    int mod_didhit[nlayers] = {0};
    int mod_shouldhit[nlayers] = {0};
    bool beam_in_range = false;
    
    while(Tscaler->GetEntry(nevent_scaler++)){
      nevent_end = -1; //Initialize end to nonsense number
    
      if(abs(beam_curr - curr_array[iQ2]) < 2){ //Check if beam curr in range
	//start event count if this is the first event in curr range
	if(!beam_in_range) nevent_start = beam_evnum; 
	beam_in_range = true;
      } 
      else { //Check if beam curr not in range
	// If the curr just dropped out then record this as the end event count
	if(beam_in_range){
	  Tscaler->GetEntry(nevent_scaler - 2); // Go back one event 
	  nevent_end = beam_evnum;
	  Tscaler->GetEntry(nevent_scaler - 1); // Set back to correct 
	}
	beam_in_range = false;
      }
      
      if(nevent_end > nevent_start){ //If we have a range now we look at data
	while(C->GetEntry(nevent++) && nevent < set_max_events){

	  bool didhit[nmodules] = {false};
	  bool shouldhit[nmodules] = {false};	  
	  
	  if(evnum > nevent_end){ //break if we are above the event range
	    nevent--;
	    break;
	  }
	  // Only anlyze events in range
	  if(evnum > nevent_start && evnum < nevent_end){ 
	    if(tr_n == 1 && ps_e > 0.15 && abs(tr_vz[0]) < 0.27){ // Physics cuts
	      for(int ilayer = 0; ilayer < nlayers; ilayer++)
		hist[ilayer]->Fill(nstrips[ilayer]);
	      
	      TLorentzVector Peprime(px[0],py[0],pz[0],p[0]);
      
	      TVector3 track_origin(tr_x[0],tr_y[0],0.0);
	      TVector3 track_dir(tr_xp[0],tr_yp[0],1.0);
      
	      for(int ihit = 0; ihit < hit_n; ihit++)
		didhit[(int)hit_mod[ihit]] = true;
      
      
	      for(int imod = 0; imod < nmodules; imod++){
		TVector3 Intersect = TrackIntersect(imod,track_origin,track_dir,iQ2);
		bool isinactivearea = IsInActiveArea(imod,Intersect,iQ2);
	
		if(isinactivearea) shouldhit[imod] = true;

		int layer;
		if(imod >= 4) layer = 4;
		else layer = imod;

		if(didhit[imod]) mod_didhit[layer]++;
		if(shouldhit[imod]) mod_shouldhit[layer]++;
	      } 
	    
	    } //End if statement on good physics cuts
	  } //End if statement on good beam current
	} //End loop over data tree
      } //End if statement on good beam current
    } //End loop over scaler tree

    for(int ilayer = 0; ilayer < nlayers; ilayer++){
      TF1 *fit = new TF1("fit","gaus");
	  
      hist[ilayer]->Fit(fit,"q0");
      
      U_occu[iQ2][ilayer] = fit->GetParameter(1) / NStripsU[ilayer];
      U_err[iQ2][ilayer] = fit->GetParameter(2) / NStripsU[ilayer];
      
      hist[ilayer]->Delete();

      GEM_eff[iQ2][ilayer] = mod_didhit[ilayer]*1.0/mod_shouldhit[ilayer];
    }
  } //End loop over Q2 runs
  

  double Q2_eff[nQ2];

  TGraphErrors *g_occu[nQ2];
  TGraph *g_eff[nQ2];
  
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TCanvas *c2 = new TCanvas("c2","",800,600);
  TLegend *legend = new TLegend(0.60,0.70,0.9,0.9);
  int icolor = 0;

  double layer_list[nlayers] = {1,2,3,4,5};

  for(int iQ2 = 0; iQ2 < nQ2; iQ2++){
  
    Q2_eff[iQ2] = probability(3,iQ2) + probability(4,iQ2) + probability(5,iQ2);
    cout<<iQ2<<" "<<Q2_eff[iQ2]<<endl;

    icolor++;

    g_eff[iQ2] = new TGraph(nlayers,layer_list,GEM_eff[iQ2]);
    g_eff[iQ2]->SetTitle("GEN GEM Layer Efficiency;Layer # ;Efficiency");
    g_eff[iQ2]->SetMarkerStyle(8);
    g_eff[iQ2]->SetMarkerColor(icolor);
    g_eff[iQ2]->SetLineColor(icolor);
    legend->AddEntry(g_eff[iQ2],Form("Q^{2} = %g",Q2_list[iQ2]),"p");

    g_occu[iQ2] = new TGraphErrors(nlayers,layer_list,U_occu[iQ2],0,U_err[iQ2]);
    g_occu[iQ2]->SetTitle("GEN GEM Layer Occupancy;Layer # ;Occupancy");
    g_occu[iQ2]->SetMarkerStyle(8);
    g_occu[iQ2]->SetMarkerColor(icolor);
    g_occu[iQ2]->SetLineColor(icolor);

    if(iQ2 == 0){
      c1->cd();
      g_eff[iQ2]->Draw("AP");
      legend->Draw("same");

      c2->cd();
      g_occu[iQ2]->Draw("AP");
      legend->Draw("same");
    }
    else{
      c1->cd();
      g_eff[iQ2]->Draw("P");

      c2->cd();
      g_occu[iQ2]->Draw("P");
    }    

    g_occu[iQ2]->GetYaxis()->SetRangeUser(0,0.2);
    g_eff[iQ2]->GetYaxis()->SetRangeUser(0.4,1.1);
    
  }


}
