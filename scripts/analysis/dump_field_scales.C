// dump_field_scales.C
// Usage: root -l -q 'dump_field_scales.C("g4sbs_output.root")'
#include <TFile.h>
#include <TObject.h>
#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TError.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

// Read a whole text file into a std::string
static std::string ReadAll(const char* path) {
  std::ifstream in(path, std::ios::in | std::ios::binary);
  if (!in) return {};
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// Extract the *last* occurrence of a numeric value following a given G4SBS UI command
// Example matches lines like: "/g4sbs/scalebbfield 0.905"
static double ExtractLastDouble(const TString& text, const char* cmd) {
  TString re = Form("%s[[:space:]]+(-?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)", cmd);
  TRegexp rx(re, /*ignoreCase*/kFALSE);

  double last = std::numeric_limits<double>::quiet_NaN();

  Ssiz_t from = 0;
  while (true) {
    Ssiz_t matchLen = 0;
    Ssiz_t pos = text.Index(rx, &matchLen, from); // pos=start, matchLen=length
    if (pos == kNPOS) break;

    TString sub = text(pos, matchLen);
    TObjArray* toks = sub.Tokenize(" \t\r\n");
    if (toks && toks->GetLast() >= 0) {
      TString val = ((TObjString*)toks->At(toks->GetLast()))->GetString();
      last = atof(val.Data());
    }
    if (toks) { toks->Delete(); delete toks; }

    from = pos + matchLen; // continue searching after this match
  }
  return last;
}

void dump_field_scales(const char* fname="g4sbs.root") {
  // Open file
  TFile f(fname, "READ");
  if (f.IsZombie()) { std::cerr << "Cannot open " << fname << "\n"; return; }

  // Grab run_data (G4SBSRunData as a TObject)
  TObject* rd = f.Get("run_data");
  if (!rd) {
    std::cerr << "No 'run_data' object found in " << fname << "\n";
    return;
  }

  // Capture rd->Print() to a temporary file
  TString tmp;
  gSystem->TempFileName(tmp);         // <-- needs a TString&
  if (tmp.IsNull()) {
    std::cerr << "Failed to create a temporary file to capture output.\n";
    return;
  }

  gSystem->RedirectOutput(tmp.Data(), "w");
  rd->Print();
  gSystem->RedirectOutput(0);         // restore stdout

  // Read captured text
  std::string meta_str = ReadAll(tmp.Data());
  gSystem->Unlink(tmp.Data());        // cleanup

  if (meta_str.empty()) {
    std::cerr << "Failed to capture run_data->Print() output.\n";
    return;
  }

  TString meta = meta_str.c_str();

  // Extract scale factors (take the *last* value if set multiple times)
  double bb_scale  = ExtractLastDouble(meta, "/g4sbs/scalebbfield");
  double sbs_scale = ExtractLastDouble(meta, "/g4sbs/scalesbsfield");

  // Report
  std::cout << "---- Field scale factors from " << fname << " ----\n";
  if (std::isnan(bb_scale))
    std::cout << "BigBite field scale (scalebbfield): not found\n";
  else
    std::cout << "BigBite field scale (scalebbfield): " << bb_scale << "\n";

  if (std::isnan(sbs_scale))
    std::cout << "SBS field scale (scalesbsfield):    not found\n";
  else
    std::cout << "SBS field scale (scalesbsfield):    " << sbs_scale << "\n";

  if (std::isnan(bb_scale) || std::isnan(sbs_scale)) {
    std::cout << "\nNotes:\n"
                 " • If 'not found', check your macros used by g4sbs include\n"
                 "   '/g4sbs/scalebbfield <val>' and '/g4sbs/scalesbsfield <val>'.\n"
                 " • The run_data Print() output contains the executed pre/post-init macros.\n";
  }
}

