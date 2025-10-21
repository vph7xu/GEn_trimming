// ────────────────────────────────────────────────────────────────
// Yield_sim.cpp     –   compute absolute MC yield per Coulomb
// gSystem->Load("libg4sbsroot.so");            // one‑time per ROOT session
// .L Yield_sim.cpp+
// double y = yield_per_C("g4sbs_LH2.root","yourCut");
// ────────────────────────────────────────────────────────────────
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include <iostream>
#include <cstring>

// g4sbs headers (make sure $G4SBS/include is on the include path)
#include "G4SBSRunData.hh"      // :contentReference[oaicite:0]{index=0}

// Yield_perC.C ------------------------------------------------------------
double Yield_perC(const char *fname, const char *cut = "")
{
  // open file + get metadata ---------------------------------------------
  std::unique_ptr<TFile> f( TFile::Open(fname,"READ") );
  G4SBSRunData *rd = nullptr;
  f->GetObject("run_data", rd);

  // sanity check
  if(!rd){ std::cerr << "run_data missing in " << fname << '\n'; return -1; }

  // pull what we need -----------------------------------------------------
  const double beamCur_uA = rd->fBeamCur;             // μA
  const double norm       = rd->fNormalization;       // events s⁻¹ per gen evt

  const double w_evt_C = norm / ( beamCur_uA * 1e-6 );// events per Coulomb

  // how many events survive your analysis cut?
  TTree *T = static_cast<TTree*>( f->Get("T") );
  const Long64_t Nacc = ( cut && *cut ) ? T->GetEntries(cut)
                                        : T->GetEntries();

  return Nacc * w_evt_C;   // absolute yield per Coulomb
}

