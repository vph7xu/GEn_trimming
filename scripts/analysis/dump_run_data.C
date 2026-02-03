// dump_run_data.C
// Usage:
//   root -l -q 'dump_run_data.C("g4sbs_output.root")'
//   root -l -q 'dump_run_data.C("g4sbs_output.root","run_data_dump.txt")'
#include <TFile.h>
#include <TObject.h>
#include <TSystem.h>
#include <TString.h>
#include <TError.h>
#include <iostream>

void dump_run_data(const char* fname="g4sbs.root", const char* outfile="") {
  // Open the file
  TFile f(fname, "READ");
  if (f.IsZombie()) { std::cerr << "Cannot open " << fname << "\n"; return; }

  // Grab run_data
  TObject* rd = f.Get("run_data");
  if (!rd) {
    std::cerr << "No 'run_data' object found in " << fname << "\n";
    f.Close();
    return;
  }

  // Optional: show the dynamic type
  std::cout << "run_data class: " << rd->IsA()->GetName() << "\n\n";

  // If an output path is provided, redirect there; otherwise print to stdout
  TString out = outfile;
  if (!out.IsNull()) {
    // Create/overwrite output file and redirect stdout
    gSystem->RedirectOutput(out.Data(), "w");
    rd->Print();        // dump everything the object chooses to print
    gSystem->RedirectOutput(0); // restore stdout
    std::cout << "Wrote run_data Print() to: " << out << "\n";
  } else {
    // Print directly to terminal
    rd->Print();        // ROOT objects typically honor Print() verbosity internally
  }

  // Optionally, also show ROOT’s generic introspection (less detailed, but sometimes useful)
  // Uncomment if you want it:
  // std::cout << "\n--- run_data->Dump() (generic ROOT dump) ---\n";
  // rd->Dump();
}

