// SumBranch.C
// Usage:
//   root -l -q 'SumBranch.C("g4sbs_elastic_*.root","T","ev.sigma","yourCut",200,"sum_hist.pdf","sum_report.txt")'

#include <TChain.h>
#include <TTreeFormula.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cmath>

void SumBranch(const char* files,
               const char* treename,
               const char* expr,
               const char* cut = "1",
               int nbins = 200,
               const char* outpdf = "SumBranchLoopHist.pdf",
               const char* outtxt = "SumBranch_report.txt",
               Long64_t maxEntries = -1)
{
  // --- open text report file ---
  std::ofstream fout(outtxt);
  if(!fout){
    std::cerr << "[ERROR] Could not open output text file: " << outtxt << std::endl;
    return;
  }

  // Helper: print to both terminal and file
  auto Tee = [&](const std::string& s){
    std::cout << s;
    fout      << s;
  };

  TChain ch(treename);
  Long64_t nfiles = ch.Add(files);
  if(nfiles <= 0){
    Tee(std::string("[ERROR] No files matched: ") + files + "\n");
    return;
  }

  Long64_t n = ch.GetEntries();
  if(maxEntries > 0 && maxEntries < n) n = maxEntries;

  TTreeFormula fCut("fCut", cut, &ch);
  TTreeFormula fExpr("fExpr", expr, &ch);

  double sum = 0.0;
  Long64_t nSel = 0;

  double xmin =  std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();

  std::vector<double> vals;
  vals.reserve(100000);

  for(Long64_t i=0; i<n; ++i){
    ch.GetEntry(i);

    // For TChain: update formulas when the underlying tree changes
    if(fCut.GetTree() != ch.GetTree()){
      fCut.UpdateFormulaLeaves();
      fExpr.UpdateFormulaLeaves();
    }

    if(fCut.EvalInstance()){
      double x = fExpr.EvalInstance();
      if(!std::isfinite(x)) continue;

      sum += x;
      nSel++;

      vals.push_back(x);
      if(x < xmin) xmin = x;
      if(x > xmax) xmax = x;
    }
  }

  std::ios::fmtflags f = std::cout.flags();
  std::streamsize    p = std::cout.precision();
  std::cout.setf(std::ios::scientific); std::cout.precision(12);
  fout.setf(std::ios::scientific);      fout.precision(12);

  Tee("\n=== SumBranchLoopHist (entry-by-entry) ===\n");
  Tee(std::string("Files matched : ") + std::to_string(nfiles) + "\n");
  Tee(std::string("Tree          : ") + treename + "\n");
  Tee(std::string("Expr          : ") + expr + "\n");
  Tee(std::string("Cut           : ") + cut + "\n");
  Tee(std::string("Entries used  : ") + std::to_string(n) + "\n");
  Tee(std::string("N selected    : ") + std::to_string(nSel) + "\n");

  {
    std::ostringstream os;
    os.setf(std::ios::scientific);
    os.precision(12);
    os << "SUM           : " << sum << "\n";
    Tee(os.str());
  }

  if(nSel <= 0){
    Tee("[WARN] No selected entries; histogram not created.\n");
    Tee(std::string("Report saved  : ") + outtxt + "\n\n");
    return;
  }

  // Handle degenerate range
  double xmin_print = xmin, xmax_print = xmax;
  if(xmin_print == xmax_print){
    double eps = (xmin_print == 0.0 ? 1.0 : std::abs(xmin_print)*0.01);
    xmin_print -= eps; xmax_print += eps;
  } else {
    double pad = 0.01*(xmax_print - xmin_print);
    xmin_print -= pad; xmax_print += pad;
  }

  {
    std::ostringstream os;
    os.setf(std::ios::scientific);
    os.precision(12);
    os << "Min(selected) : " << xmin_print << "\n";
    os << "Max(selected) : " << xmax_print << "\n";
    Tee(os.str());
  }

  TString hname = TString::Format("h_%s",
      TString(expr).ReplaceAll(".","_").ReplaceAll("[","_").ReplaceAll("]","_").Data());

  TH1D *h = new TH1D(hname, TString::Format("%s; %s; Counts", expr, expr),
                     nbins, xmin_print, xmax_print);

  for(double x : vals) h->Fill(x);

  TCanvas *c = new TCanvas("c_sumbranch", "SumBranchLoopHist", 900, 650);
  h->SetLineWidth(2);
  h->Draw("hist");
  c->SaveAs(outpdf);

  Tee(std::string("Histogram saved: ") + outpdf + "\n");
  Tee(std::string("Report saved   : ") + outtxt + "\n\n");

  // restore cout formatting
  std::cout.flags(f);
  std::cout.precision(p);

  fout.close();
}

