//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   This script will take all the optics and H2 data and plot
//   the important correlations.
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

//Define limits for different plots
double Beammin = -3, Beammax = 4;
double rastermin = 28000, rastermax = 60000;
double pmin = -10, pmax = 10;
double dxmin_GEN3 = -2.2, dxmax_GEN3 = -1.5;
double dxmin_GEN2 = -3.5, dxmax_GEN2 = -2;
double dymin = -0.4, dymax = 0.4;
double ytarmin = -30, ytarmax = 30;
double xptarmin = -20, xptarmax = 20;
double yptarmin = -20, yptarmax = 20;
double ztarmin = -30, ztarmax = 30;
double W2max = 1.2;


struct data_info {
  TString data;
  TString config;
  TString var;
  TString optics;
  TString title;
  TString units;
  double ymin;
  double ymax;
};


data_info SetData(TString input){

  data_info info;

  //Input is expected to be like this, "H2 dx GEN2 new"
  //The last value is the optics matrix type
  istringstream iss(input.Data());
  string token;
  getline(iss,token,' ');
  TString data = token;
  getline(iss,token,' ');
  TString var = token;
  getline(iss,token,' ');
  TString config = token;
  getline(iss,token,' ');
  TString optics = token;

  if(data != "optics" && data != "H2"){
    cout<<"Data must be type 'optics' or 'H2'"<<endl;
    exit(1);
  }


  info.data = data;
  info.config = config;
  info.var = "(trP/pcentral - 1)*100";
  info.optics = optics;
  info.title = "#deltap";
  info.units = " (%)";
  info.ymin = pmin;
  info.ymax = pmax;

  //Set the variables for this data set
  if(var == "dx"){
    info.var = "dx";
    info.title = "#Deltax";
    info.units = " (m)";
    info.ymin = dxmin_GEN2;
    info.ymax = dxmax_GEN2;
    if(config == "GEN3"){
      info.ymin = dxmin_GEN3;
      info.ymax = dxmax_GEN3;
    }
  }
  if(var == "dy"){
    info.var = "dy";
    info.title = "#Deltay";
    info.units = " (m)";
    info.ymin = dymin;
    info.ymax = dymax;
  }
  if(var == "ytar"){
    info.var = "ytar";
    info.title = "#Deltay_{tg}";
    info.units = " (mm)";
    info.ymin = ytarmin;
    info.ymax = ytarmax;
  }
  if(var == "xptar"){
    info.var = "xptar";
    info.title = "#Delta#theta_{tg}";
    info.units = " (mrad)";
    info.ymin = xptarmin;
    info.ymax = xptarmax;
  }
  if(var == "yptar"){
    info.var = "yptar";
    info.title = "#Delta#phi_{tg}";
    info.units = " (mrad)";
    info.ymin = yptarmin;
    info.ymax = yptarmax;
  }
  if(var == "ztar"){
    info.var = "ztar";
    info.title = "#Deltaz";
    info.units = " (mm)";
    info.ymin = ztarmin;
    info.ymax = ztarmax;
  }
  
  return info;
}
	     
////// Function for fitting 1D histograms with a gaussian ///////////////////
void draw_gaus(TH1D *hist, TLegend *leg){

  TF1 *fit1 = new TF1("g1","gaus");
  hist->Fit("g1","qN");

  TF1 *fit2 = new TF1("g2","gaus",fit1->GetParameter(1) - 1.5*fit1->GetParameter(2),fit1->GetParameter(1) + 1.5*fit1->GetParameter(2));
  hist->Fit("g2","qNR");
  
  leg->AddEntry((TObject*)0,Form("Gaus Fit #sigma = %.3g",fit2->GetParameter(2)),"");
}


// Function used to fit a linear correlation on a plot
void draw_corr(TH2F *hist_orig, TH2F *hist_new = NULL, vector<double> lincuts = {}){

  TH2F *hist = hist_orig;
  TH2F *hist_yrange = (TH2F*)hist->Clone("temp");
  
  //Get axis values  
  double xmin = hist->GetXaxis()->GetXmin();
  double xmax = hist->GetXaxis()->GetXmax();
  
  hist_yrange->Reset();

  //Only do this if linear cuts were defined
  if(lincuts.size() != 0){

    double slope_low = (lincuts[1] - lincuts[0])/(xmax - xmin);
    double slope_high = (lincuts[3] - lincuts[2])/(xmax - xmin);
    double b_low = lincuts[1] - slope_low*xmax;
    double b_high = lincuts[3] - slope_high*xmax;

    for (Int_t binX = 1; binX <= hist->GetNbinsX(); binX++)
      {
	for (Int_t binY = 1; binY <= hist->GetNbinsY(); binY++)
	  {
	    double y = hist->GetYaxis()->GetBinCenter(binY);
	    double x = hist->GetXaxis()->GetBinCenter(binX);
	    //Fill a new histogram only if it is within the linear cuts
	    if(y > slope_low*x + b_low && y < slope_high*x + b_high)
	      hist_yrange->SetBinContent(binX,binY,hist->GetBinContent(binX,binY));
	    
	  }
      }
  
    TLine *lin_low = new TLine(xmin,lincuts[0],xmax,lincuts[1]);
    lin_low->SetLineColor(kRed);  
    //lin_low->Draw("same");

    TLine *lin_high = new TLine(xmin,lincuts[2],xmax,lincuts[3]);
    lin_high->SetLineColor(kRed);  
    //lin_high->Draw("same");
    
    //Set this equal tothe new histogram with the linear cuts
    hist = (TH2F*)hist_yrange->Clone("temp2");
  } //end on linear cuts
  
  TF1 *linfit = new TF1("f1","pol1",xmin, xmax);
  //Fit a linear line
  hist->Fit("f1","qN");
  linfit->Draw("same");

  //Loop over the original histogram
  for (Int_t binX = 1; binX <= hist_orig->GetNbinsX(); binX++)
    {
      for (Int_t binY = 1; binY <= hist_orig->GetNbinsY(); binY++)
        {
	  double x = hist_orig->GetXaxis()->GetBinCenter(binX);
	  double y = hist_orig->GetYaxis()->GetBinCenter(binY);

	  double binContent = hist_orig->GetBinContent(binX, binY);

	  //If this histogram was defined then we use it to remove the linear correlation
	  if(hist_new != NULL){
	    double y_new = y - (linfit->GetParameter(1)*x + linfit->GetParameter(0));
	    int binY_new = hist_new->GetYaxis()->FindBin(y_new);
	    hist_new->SetBinContent(binX,binY_new,binContent);
	  }
	}
    }
  
  double corr_coeff = hist_orig->GetCorrelationFactor();

  //Write the correlation coefficient on the plot
  TPaveText *pt = new TPaveText(.15,.15,.4,.25,"ndc");
  pt->AddText(Form("Corr Coeff = %.2g",corr_coeff));
  pt->SetFillColor(0);
  pt->Draw("same");

}

//Function used to plot the beam variables, BPMs and rasters
void Beam_variables(){

  TString filename = "../outfiles/QE_data_GEN3_sbs100p_nucleon_p_model1.root";
  TFile *file = new TFile(filename,"read");
  TTree *T = (TTree*)file->Get("Tout");

  TH2F *hhBPM = new TH2F("hhBPM","BPM Projections to Target;Beam x (mm);Beam y (mm)",100,Beammin,Beammax,100,Beammin,Beammax);
  TH2F *hhRaster = new TH2F("hhRaster","Raster Current Y vs X;Raster Current X;Raster Current Y",100,rastermin,rastermax,100,rastermin,rastermax);
  TH2F *hhBeam = new TH2F("hhBeam","Beam Y vs X;Beam Pos X (mm);Beam Pos Y (mm)",100,Beammin,Beammax,100,Beammin,Beammax);
  TH2F *hhBPMx_Rasterx = new TH2F("hhBPMx_Rasterx","BPMA X vs Raster X;Raster Current X;BPMA X (mm)",100,rastermin,rastermax,100,Beammin,Beammax);
  TH2F *hhBPMy_Rastery = new TH2F("hhBPMy_Rastery","BPMA Y vs Raster Y;Raster Current Y;BPMA Y (mm)",100,rastermin,rastermax,100,Beammin,Beammax);

  
  TCut runcut = "runnum < 2800";

  T->Draw("BPMy:BPMx>>hhBPM",runcut,"goff");
  T->Draw("Rastery:Rasterx>>hhRaster",runcut,"goff");
  T->Draw("vy*1000:vx*1000>>hhBeam",runcut,"goff");
  T->Draw("BPMAx*1000:Rasterx>>hhBPMx_Rasterx",runcut,"goff");
  T->Draw("BPMAy*1000:Rastery>>hhBPMy_Rastery",runcut,"goff");


  //////// Beam Position Stability over multiple runs /////////////////////
  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->Divide(2,3);

  c1->cd(1);
  T->Draw("BPMAx:runnum","","colz");
  TH1F *h = (TH1F*)gPad->GetPrimitive("htemp");
  h->SetTitle("BPMA X vs Run Number;Run Number;BPMA X");
  
  c1->cd(2);
  T->Draw("BPMAy:runnum","","colz");
  TH1F *h2 = (TH1F*)gPad->GetPrimitive("htemp");
  h2->SetTitle("BPMA Y vs Run Number;Run Number;BPMA Y"); 
 
  c1->cd(3);
  T->Draw("Rasterx:runnum","","colz");
  TH1F *h3 = (TH1F*)gPad->GetPrimitive("htemp");
  h3->SetTitle("Raster Current X vs Run Number;Run Number;Raster Current X"); 

  c1->cd(4);
  T->Draw("Rasterx:runnum","","colz");
  TH1F *h4 = (TH1F*)gPad->GetPrimitive("htemp");
  h4->SetTitle("Raster Current Y vs Run Number;Run Number;Raster Current Y"); 

  c1->cd(5);
  T->Draw("vx*1000:runnum","","colz");
  TH1F *h5 = (TH1F*)gPad->GetPrimitive("htemp");
  h5->SetTitle("Beam Position X vs Run Number;Run Number;Beam Position X (mm)"); 

  c1->cd(6);
  T->Draw("vy*1000:runnum","","colz");
  TH1F *h6 = (TH1F*)gPad->GetPrimitive("htemp");
  h6->SetTitle("Beam Position Y vs Run Number;Run Number;Beam Position Y (mm)"); 

  ////////////////////////////////////////////////////////////////////

  //////// Beam X/Y plots for BPM/Raster/Vertex /////////////////////
  TCanvas *c2 = new TCanvas("c2","",800,800);
  c2->Divide(2,3);

  TPaveText *pt = new TPaveText(.45,.72,.88,.88,"ndc");
  pt->AddText("Cut on runs with similar beam pos");
  pt->AddText("for all future plots");
  pt->SetFillColor(0);

  c2->cd(1);
  hhBPM->Draw("colz");
  pt->Draw("same");

  c2->cd(2);
  hhRaster->Draw("colz");
  
  c2->cd(3);
  hhBeam->Draw("colz");

  c2->cd(5);
  hhBPMx_Rasterx->Draw("colz");

  c2->cd(6);
  hhBPMy_Rastery->Draw("colz");

  //////////////////////////////////////////////////////////////////////
  
  TString outfile = "../../plots/Beam_Variables.pdf";
  cout<<"PDF: "<<outfile<<endl;  

  c1->Print(outfile + "(");
  c2->Print(outfile + ")");

}

//This function will automatically draw a histogram with the beam x/y correlation for the variable input
void make_canvas(TString input, TTree *T,TCanvas *c,TH1D *hist[] = {}){

  data_info info = SetData(input); //Get the info for this data type
  
  //Define beam x/y histograms
  TH2F *hhBeamx = new TH2F("hhBeamx",info.title + " vs Beam X correlation;Beam x (mm);" + info.title + info.units,100,Beammin,Beammax,100,info.ymin,info.ymax);
  TH2F *hhBeamy = new TH2F("hhBeamy",info.title + " vs Beam Y correlation;Beam y (mm);" + info.title + info.units,100,Beammin,Beammax,100,info.ymin,info.ymax);
  TH2F *hhBeamx_new = new TH2F("hhBeamx_new",info.title + " vs Beam X correlation;Beam x (mm);" + info.title + info.units,100,Beammin,Beammax,100,info.ymin,info.ymax);
  TH2F *hhBeamy_new = new TH2F("hhBeamy_new",info.title + " vs Beam Y correlation;Beam y (mm);" + info.title + info.units,100,Beammin,Beammax,100,info.ymin,info.ymax);

  TCut total_cut = "";
  
  //If H2 data we must add some extra cuts
  if(info.data == "H2"){
    TCut run_cut = "runnum > 2008";     //Cutting out some runs that seem to have missed the target for GEN2
    TCut W2cut = Form("W2 < %g",W2max);
    total_cut = run_cut && "pCut";
  }


  if(info.data == "optics"){
    T->Draw("(" + info.var + "-" + info.var + "T)*1000:targx*1000>>hhBeamx",total_cut,"goff");
    T->Draw("(" + info.var + "-" + info.var + "T)*1000:targy*1000>>hhBeamy",total_cut,"goff");
  }

  if(info.data == "H2"){
    T->Draw(info.var + ":vx*1000>>hhBeamx",total_cut,"goff");
    T->Draw(info.var + ":vy*1000>>hhBeamy",total_cut,"goff");
  }

  if(info.data == "H2"){
    
    TPaveText *pt = new TPaveText(.65,.72,.88,.88,"ndc");
    pt->AddText("2#sigma HCal Spot Cut");
    pt->AddText("On all plots here");
    pt->SetFillColor(0);
  
    c->cd(1);
    hhBeamx->Draw("colz");
    draw_corr(hhBeamx);
    pt->Draw("same");    

    c->cd(2);
    hhBeamy->Draw("colz");

    //For deltap data we remove the correlation and check the result
    if(info.title == "#deltap"){
      vector<double> lincuts = {1,-7,8,1};  //Number for GEN2 H2
      if(info.config == "GEN2" && info.optics == "old") draw_corr(hhBeamy,hhBeamy_new,lincuts);
      else draw_corr(hhBeamy,hhBeamy_new);
      

      hist[0] = hhBeamy->ProjectionY();
      hist[1] = hhBeamy_new->ProjectionY();
    
      hist[0]->SetName("notcorr_"+info.title);
      hist[0]->SetTitle(info.title + " Comparison");

      hist[1]->SetName("corr_"+info.title);
      hist[1]->SetTitle(info.title + " Comparison");
      hist[1]->SetLineColor(kRed);
      
    }
    else draw_corr(hhBeamy);
     
  }

  if(info.data == "optics"){

    c->cd(1);
    hhBeamx->Draw("colz");
    draw_corr(hhBeamx,hhBeamx_new);

    c->cd(2);
    hhBeamy->Draw("colz");
    draw_corr(hhBeamy,hhBeamy_new);
  
    //For ytar/yptar we remove the beamx correltation
    if(info.var == "ytar" || info.var == "yptar"){
      hist[0] = hhBeamx->ProjectionY();
      hist[1] = hhBeamx_new->ProjectionY();
      //hhBeamx_new->Draw("colz");
      //draw_corr(hhBeamx_new);    
    }
    //For xptar we remove the beamy correltation
    if(info.var == "xptar"){
      hist[0] = hhBeamy->ProjectionY();
      hist[1] = hhBeamy_new->ProjectionY();
      //hhBeamy_new->Draw("colz");
      //draw_corr(hhBeamy_new);    
    }

    hist[0]->SetName("notcorr_"+info.title);
    hist[0]->SetTitle(info.title + " Comparison");

    hist[1]->SetName("corr_"+info.title);
    hist[1]->SetTitle(info.title + " Comparison");
    hist[1]->SetLineColor(kRed);

  }

}

// This is used to make simple 1D plot comparisons between two files with the same variable
void make_comparison_1D(TTree *T_nobeam,TTree *T_beam,TString input){

  data_info info = SetData(input);

  TString settitle = info.title + " Matrix Comparison (Optics Data);" + info.title + info.units + ";";
  if(info.data == "H2") settitle = info.title + " Matrix Comparison (H2 Data);" + info.title + info.units + ";";

  TH1D *hnobeam = new TH1D("hnobeam",settitle,100,info.ymin,info.ymax);
  TH1D *hbeam = new TH1D("hbeam",settitle,100,info.ymin,info.ymax);
  TH1D *hbeam_new = new TH1D("hbeam_new",settitle,100,info.ymin,info.ymax);

  TPaveText *pt = new TPaveText(.65,.72,.88,.88,"ndc");
  pt->SetFillColor(0);

  TCut total_cut = "";
  
  if(info.data == "H2"){
    TCut run_cut = "runnum > 2008";     //Cutting out some runs that seem to have missed the target for GEN2
    TCut W2cut = Form("W2 < %g",W2max);
 
    if(info.title == "#deltap"){
      total_cut = run_cut && "pCut";
      pt->AddText("2#sigma HCal Spot Cut");
    } else {
      total_cut = run_cut && W2cut; 
      pt->AddText(Form("W^{2} < %g GeV",W2max));  
    }
    
  }
  
  //Draw the variable on the same plot and normalize the y axis
  if(info.data == "optics"){
    T_nobeam->Draw("(" + info.var + "-" + info.var + "T)*1000>>hnobeam",total_cut,"goff");
    T_beam->Draw("(" + info.var + "-" + info.var + "T)*1000>>hbeam",total_cut,"goff");

    hnobeam->Scale(1/hnobeam->Integral());
    hbeam->Scale(1/hbeam->Integral());
  }
  //Draw the variable on the same plot and normalize the y axis
  if(info.data == "H2"){
    T_nobeam->Draw(info.var + ">>hnobeam",total_cut,"goff");
    T_beam->Draw(info.var + ">>hbeam",total_cut,"goff");
    
    hnobeam->Scale(1/hnobeam->Integral());
    hbeam->Scale(1/hbeam->Integral());
   
  }

  TLegend *leg = new TLegend(0.11,0.72,0.44,0.89);
  leg->SetLineColor(0);
  
  leg->AddEntry(hnobeam,"Old Optimization","l"); 
  if(info.var == "xptar" || info.var == "yptar") leg->AddEntry((TObject*)0,Form("Hist RMS = %g",hnobeam->GetRMS()),""); 
  else draw_gaus(hnobeam,leg);
  
  leg->AddEntry(hbeam,"New Optimization","l");  
  if(info.var == "xptar" || info.var == "yptar")leg->AddEntry((TObject*)0,Form("Hist RMS = %g",hbeam->GetRMS()),"");  
  else draw_gaus(hbeam,leg);
  
  
  hbeam->Draw("hist");
  hnobeam->Draw("same hist");
  hbeam->SetLineColor(kRed);
  leg->Draw("same");
  if(info.data == "H2") pt->Draw("same");

}

//For H2 data we are interested in the deltap and dx/dy plots
void H2_correlations(TString config = "GEN3", TString optics = "old")
{
  gStyle->SetOptStat(0);
  
  TString filename = "../outfiles/QE_data_" + config + "_sbs100p_nucleon_p_model1.root";
  if(optics == "new") filename = "../outfiles/QE_new_matrix_" + config + "_sbs100p_nucleon_p_model1.root";
  TFile *file = new TFile(filename,"read");
  TTree *T = (TTree*)file->Get("Tout");
  
  TH1D *hcorrect[5][2] = {};

  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->Divide(2,1);
  //Draw 2D deltap plot and get the corrected histogram
  make_canvas("H2 trP " + config + " " + optics,T,c1,hcorrect[0]);

  TLegend *leg2 = new TLegend(0.11,0.72,0.44,0.89);
  leg2->SetLineColor(0);  
  
  TCanvas *c2 = new TCanvas("c2","",800,600);
  //Draw the corrected and uncorrected deltap 1D plot
  hcorrect[0][1]->Draw();
  hcorrect[0][0]->Draw("same");

  leg2->AddEntry(hcorrect[0][0],"With Correlation","l");  
  draw_gaus(hcorrect[0][0],leg2);
  
  leg2->AddEntry(hcorrect[0][1],"Correlation Removed","l");
  draw_gaus(hcorrect[0][1],leg2);
  leg2->Draw("same");
  
  TCanvas *c3 = new TCanvas("c3","",1200,600);
  c3->Divide(2,1);
  //Draw the dx plot
  make_canvas("H2 dx " + config + " " + optics,T,c3);
  
  TCanvas *c4 = new TCanvas("c4","",1200,600);
  c4->Divide(2,1);
  //Draw the dy plot
  make_canvas("H2 dy " + config + " " + optics,T,c4);
  
  TString outfile = "../../plots/H2_beam_correlations_" + config + "_" + optics + "_matrix.pdf";
  cout<<"PDF: "<<outfile<<endl;

  c1->Print(outfile + "(");
  c2->Print(outfile);
  c3->Print(outfile);
  c4->Print(outfile + ")");
  
}

//Takes the optics data and draws all the correlations plots
void Optics_correlations(TString kin = "GEN3",TString config = "original"){

  // fileid = -3 (original), -2 (new), -1 (old) matrix 
  int fileid = -3;
  TString info1 = "Data replayed using original optics matrix";
  TString info2 = "For this data the 'true' value does NOT include the beam position";

  TString cfg_lower = "gen3";
  
  if(kin == "GEN4") cfg_lower = "gen4";
  else if(kin == "GEN2") cfg_lower = "gen4";
  else if(kin == "GEN4b") cfg_lower = "gen4b";
  TString matrix_name = "optics_" + cfg_lower + ".txt";
  

  //Adds some description on the pdf of what this data set means
  if(config == "old"){
    fileid = -1;
    info1 = "Data replayed using original optics matrix";
    info2 = "For this data the 'true' value DOES include the beam position";
  }
  if(config == "new"){
    fileid = -2;
    matrix_name = "optics_" + cfg_lower + "_beam.txt";
    if(kin == "GEN4b") matrix_name = "optics_" + cfg_lower + "_new.txt";
   
    info1 = "Data replayed using new optics matrix, including the beam position";
    info2 = "For this data the 'true' value DOES include the beam position";
  }

  TString filename = Form("/w/halla-scshelf2102/sbs/jeffas/genOpt/hist/Optics_gen3_7foil_%i_fit_tree.root",fileid);

  if(kin == "GEN4")
    filename = Form("/w/halla-scshelf2102/sbs/jeffas/genOpt/hist/Optics_gen4_7foil_%i_fit_tree.root",fileid);
  else if(kin == "GEN2")
    filename = Form("/w/halla-scshelf2102/sbs/jeffas/genOpt/hist/Optics_gen2_7foil_%i_fit_tree.root",fileid);
  else if(kin == "GEN4b")
    filename = Form("/w/halla-scshelf2102/sbs/jeffas/genOpt/hist/Optics_gen4b_7foil_%i_fit_tree.root",fileid);
  

  TFile *file = new TFile(filename,"read");
  TTree *T = (TTree*)file->Get("TFit");

  ifstream matrix(matrix_name);

  //Histograms which will show results after correlation corrections
  TH1D *hcorrect[3][2] = {};

  //Draw all the 2D correlations and get the corrected 1D plots
  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->Divide(2,1);
  make_canvas("optics ytar",T,c1,hcorrect[0]);

  TCanvas *c2 = new TCanvas("c2","",1200,600);
  c2->Divide(2,1);
  make_canvas("optics yptar",T,c2,hcorrect[1]);

  TCanvas *c3 = new TCanvas("c3","",1200,600);
  c3->Divide(2,1);
  make_canvas("optics xptar",T,c3,hcorrect[2]);
  
  TCanvas *c4 = new TCanvas("c4","",1200,600);
  c4->Divide(2,2);

  ////////// y_tg result comparison //////////
  c4->cd(1);
  hcorrect[0][0]->Draw();  
  hcorrect[0][1]->Draw("same");
  

  TLegend *leg2 = new TLegend(0.11,0.72,0.47,0.89);
  leg2->SetLineColor(0);

  leg2->AddEntry(hcorrect[0][0],"With Correlation","l");
  draw_gaus(hcorrect[0][0],leg2);

  leg2->AddEntry(hcorrect[0][1],"Correlation Removed","l");
  draw_gaus(hcorrect[0][1],leg2); 
  leg2->Draw("same");
  ////////////////////////////////////////////////

  ////////// phi_tg result comparison //////////
  c4->cd(2);
  hcorrect[1][1]->Draw();
  hcorrect[1][0]->Draw("same");

  TLegend *leg3 = new TLegend(0.11,0.72,0.43,0.89);
  leg3->SetLineColor(0);

  leg3->AddEntry(hcorrect[0][0],"With Correlation","l");
  draw_gaus(hcorrect[1][0],leg3);

  leg3->AddEntry(hcorrect[0][1],"Correlation Removed","l");
  draw_gaus(hcorrect[1][1],leg3); 
  leg3->Draw("same");
  ////////////////////////////////////////////////

  ////////// theta_tg result comparison //////////
  c4->cd(3);
  hcorrect[2][1]->Draw();
  hcorrect[2][0]->Draw("same");

  TLegend *leg4 = new TLegend(0.11,0.72,0.43,0.89);
  leg4->SetLineColor(0);

  leg4->AddEntry(hcorrect[0][0],"With Correlation","l");
  draw_gaus(hcorrect[1][0],leg4);

  leg4->AddEntry(hcorrect[0][1],"Correlation Removed","l");
  draw_gaus(hcorrect[1][1],leg4); 
  leg4->Draw("same");

  ////////////////////////////////////////////////

  TCanvas *c5 = new TCanvas("c5","",1200,800);
  TPaveText *pt = new TPaveText(0.01,0.01,0.99,0.99,"ndc");
  pt->AddText(info1 + ". All #Delta's defined as 'reconstructed' - 'true' value");
  pt->AddText(info2);
  pt->AddText("");
  pt->AddText("");
  pt->AddText("Matrix:");

  std::string line;
  int iline = 0;
  while (std::getline(matrix, line)){
    iline++;
    if(iline == 1) continue;
    pt->AddText(line.c_str());
  }

  pt->SetFillColor(0);
  pt->Draw();
  

  TString outfile = "../../plots/optics_correlations_" + kin + "_" + config + "_matrix.pdf";
  cout<<"PDF: "<<outfile<<endl;  

  c5->Print(outfile + "(");
  c1->Print(outfile);
  c2->Print(outfile);
  c3->Print(outfile);
  c4->Print(outfile + ")");

}

// Creates a simple summary of how the resolutions compare between the old and new matrix
void Matrix_comparison(TString config = "GEN3"){

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  TCanvas *c2 = new TCanvas("c2","",800,600);
  TCanvas *c3 = new TCanvas("c3","",800,600);
  TCanvas *c4 = new TCanvas("c4","",800,600);

  //Right now we do not have optics for GEN2 for rastered beam
  //if(config == "GEN3" || config == "GEN4"){

    TString cfg_lower = "gen3";

    if(config == "GEN2") cfg_lower = "gen2";
    if(config == "GEN4") cfg_lower = "gen4";
    if(config == "GEN4b") cfg_lower = "gen4b";

    //Recall that -1 = old matrix and -2 = new matrix
    TFile *file_nobeam = new TFile("/w/halla-scshelf2102/sbs/jeffas/genOpt/hist/Optics_" + cfg_lower + "_7foil_-1_fit_tree.root","read");
    TFile *file_beam = new TFile("/w/halla-scshelf2102/sbs/jeffas/genOpt/hist/Optics_" + cfg_lower + "_7foil_-2_fit_tree.root","read");
    TTree *T_nobeam = (TTree*)file_nobeam->Get("TFit");
    TTree *T_beam = (TTree*)file_beam->Get("TFit");
    
    ///////// Compare optics run reconstruction //////////////////////
    c1->Divide(2,2);

    c1->cd(1);
    make_comparison_1D(T_nobeam,T_beam,"optics ytar " + config);

    c1->cd(2);
    make_comparison_1D(T_nobeam,T_beam,"optics yptar " + config);

    c1->cd(3);
    make_comparison_1D(T_nobeam,T_beam,"optics xptar " + config);

    c1->cd(4);
    make_comparison_1D(T_nobeam,T_beam,"optics ztar " + config);

    TH2F *hhsieve_xy = new TH2F("hhsieve_xy","Sieve Reconstruction;Sieve Y (m);Sieve X (m)",100,-0.15,0.15,100,-0.35,0.35);
    TH1F *hsieve_vz_old = new TH1F("hsieve_vz_old","Optics Foils Vertex;z (m)",200,-0.3,0.3);
    TH1F *hsieve_vz_new = new TH1F("hsieve_vz_new","Optics Foils Vertex;z (m)",200,-0.3,0.3);

    c3->cd();
    T_beam->Draw("xs:ys>>hhsieve_xy","","colz");

    c4->cd();
    T_beam->Draw("ztar>>hsieve_vz_new");
    T_nobeam->Draw("ztar>>hsieve_vz_old","","same");
    hsieve_vz_new->SetLineColor(kRed);

    TLegend *leg = new TLegend(0.60,0.72,0.89,0.89);
    leg->SetLineColor(0);
    leg->AddEntry(hsieve_vz_old,"Old Optimization","l"); 
    leg->AddEntry(hsieve_vz_new,"New Optimization","l");  
    leg->Draw("same");

    //}
  ////////////////////////////////////////////////////////////////
  
 
    ///////// Compare H2 Runs with different optics ////////////////
    TString H2_filename_nobeam = "../outfiles/QE_old_matrix_" + config + "_sbs100p_nucleon_p_model1.root";
    TString H2_filename_beam = "../outfiles/QE_new_matrix_" + config + "_sbs100p_nucleon_p_model1.root";
    TFile *H2_file_nobeam = new TFile(H2_filename_nobeam,"read");
    TFile *H2_file_beam = new TFile(H2_filename_beam,"read");
    TTree *T_H2_nobeam = (TTree*)H2_file_nobeam->Get("Tout");
    TTree *T_H2_beam = (TTree*)H2_file_beam->Get("Tout");  
 
    c2->cd();
    make_comparison_1D(T_H2_nobeam,T_H2_beam,"H2 trP " + config);

 
   ////////////////////////////////////////////////////////////////
  
  TString outfile = "../../plots/Matrix_comparison_" + config + ".pdf";
  cout<<"PDF: "<<outfile<<endl;

  if(config == "GEN3" || config == "GEN4b"){
    c1->Print(outfile + "(");
    c2->Print(outfile);
    c3->Print(outfile);
    c4->Print(outfile + ")");
  } 
  else if(config == "GEN4"){
    c1->Print(outfile + "(");
    c3->Print(outfile);
    c4->Print(outfile + ")");
  }
  else if(config == "GEN2"){
    c2->Print(outfile);
  }

}

//This is the main function
void beam_correlations(){

  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kError;

  //Beam_variables();
  /*
  H2_correlations("GEN2","old");
  H2_correlations("GEN2","new");
  H2_correlations("GEN3","old"); 
  H2_correlations("GEN3","new");
  Optics_correlations("GEN3","old");
  Optics_correlations("GEN3","new");
  Optics_correlations("GEN4","old");
  Optics_correlations("GEN4","new");
  Matrix_comparison("GEN2");
  Matrix_comparison("GEN3");
  Matrix_comparison("GEN4");
  */

  //Optics_correlations("GEN4b","old");
  //Optics_correlations("GEN4b","new");
  Matrix_comparison("GEN4b");
}
