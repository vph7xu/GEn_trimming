#!/usr/bin/env python3
"""Patch GEn_trimming/scripts/replay/QuasiElastic_ana_withsbsgems.C
for dual SBSOFF/SBSON data replay input.

The patched macro:
  * uses SBSOFF as the primary event chain (BB/HCAL/kinematics),
  * uses SBSON only as a per-run companion chain for SBS GEM/SBS tracks,
  * matches SBSOFF and SBSON by CODA event number within each run,
  * never loads simulation through LoadRawRootFiles(..., 0),
  * keeps events with no SBS track,
  * records run_has_sbson, sbson_event_match, and has_sbs_track flags,
  * writes NaN for SBS-track quantities when no track is available.

Usage:
    python3 apply_dual_replay_update_withsbsgems_v3.py /path/to/GEn_trimming
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


def replace_once(text: str, old: str, new: str, label: str) -> str:
    """Replace one block, tolerating trailing spaces/tabs on source lines."""
    count = text.count(old)
    if count == 1:
        return text.replace(old, new, 1)
    if count > 1:
        raise RuntimeError(
            f"Expected exactly one occurrence for {label}, found {count}."
        )

    # Fallback: the repository has had whitespace-only differences on blank
    # lines. Match line-for-line after removing only trailing spaces/tabs,
    # while leaving the rest of the macro untouched.
    src_lines = text.splitlines(keepends=True)
    old_lines = old.splitlines()
    old_norm = [line.rstrip(" \t") for line in old_lines]

    matches = []
    n = len(old_lines)
    for i in range(0, len(src_lines) - n + 1):
        candidate = [
            line.rstrip("\r\n").rstrip(" \t")
            for line in src_lines[i:i+n]
        ]
        if candidate == old_norm:
            matches.append(i)

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one occurrence for {label}, found {len(matches)} "
            "even after ignoring trailing whitespace. Your local macro likely "
            "contains substantive changes in that section."
        )

    i = matches[0]
    replacement = new.splitlines(keepends=True)
    src_lines[i:i+n] = replacement
    return "".join(src_lines)


def main() -> int:
    if len(sys.argv) != 2:
        print(f"Usage: {Path(sys.argv[0]).name} /path/to/GEn_trimming", file=sys.stderr)
        return 2

    repo = Path(sys.argv[1]).expanduser().resolve()
    macro = repo / "scripts" / "replay" / "QuasiElastic_ana_withsbsgems.C"
    if not macro.is_file():
        raise FileNotFoundError(f"Could not find {macro}")

    original = macro.read_text(encoding="utf-8")
    text = original

    # ------------------------------------------------------------------
    # Includes
    # ------------------------------------------------------------------
    text = replace_once(
        text,
        '#include <iostream>\n',
        '#include <iostream>\n#include <algorithm>\n#include <cmath>\n#include <limits>\n',
        "standard-library includes",
    )

    # ------------------------------------------------------------------
    # Function signature: optional explicit SBSOFF/SBSON directories.
    # The normal config file should still be the non-SBS-tracking-biased one.
    # ------------------------------------------------------------------
    old_sig = (
        'int QuasiElastic_ana_withsbsgems(const std::string configfilename, '
        'std::string filebase="/volatile/halla/sbs/vimukthi/outfiles/pass3/H2/SBSON/QE_data")\n'
    )
    new_sig = (
        'int QuasiElastic_ana_withsbsgems(\n'
        '    const std::string configfilename,\n'
        '    std::string filebase="/volatile/halla/sbs/vimukthi/outfiles/pass3/H2/SBSON/QE_data",\n'
        '    std::string sbsoff_dir="",\n'
        '    std::string sbson_dir="")\n'
    )
    text = replace_once(text, old_sig, new_sig, "function signature")

    # ------------------------------------------------------------------
    # Resolve the two DATA directories before building the primary chains.
    # ------------------------------------------------------------------
    old_load = '''  Utilities::KinConf kin_info = Utilities::LoadKinConfig(configdir + configfilename,1);\n\n  getDB(kin_info.conf);\n\n  // parsing trees\n  TChain *C = LoadRawRootFiles(kin_info, 1);\n  TChain *E = LoadRawRootFiles_E(kin_info, 1);\n'''

    new_load = '''  Utilities::KinConf kin_info = Utilities::LoadKinConfig(configdir + configfilename,1);\n\n  // ------------------------------------------------------------------\n  // Dual-replay DATA input\n  //   SBSOFF = authoritative event sample for BB/HCAL/kinematics\n  //   SBSON  = SBS GEM/SBS tracking companion for runs that exist there\n  // ------------------------------------------------------------------\n  auto replace_dir_token = [](std::string path, const std::string &from,\n                              const std::string &to) {\n    const std::size_t pos = path.find(from);\n    if (pos != std::string::npos) path.replace(pos, from.size(), to);\n    return path;\n  };\n\n  if (sbsoff_dir.empty()) {\n    sbsoff_dir = kin_info.rootfile_dir;\n    if (sbsoff_dir.find("SBSON") != std::string::npos)\n      sbsoff_dir = replace_dir_token(sbsoff_dir, "SBSON", "SBSOFF");\n  }\n\n  if (sbson_dir.empty()) {\n    sbson_dir = sbsoff_dir;\n    if (sbson_dir.find("SBSOFF") != std::string::npos)\n      sbson_dir = replace_dir_token(sbson_dir, "SBSOFF", "SBSON");\n  }\n\n  if (sbsoff_dir == sbson_dir) {\n    std::cerr << "ERROR: could not determine distinct SBSOFF and SBSON directories.\\n"\n              << "Pass them explicitly as the 3rd and 4th macro arguments." << std::endl;\n    return 2;\n  }\n\n  // An SBS-track requirement in the global cut would bias a tracking-efficiency\n  // denominator and would also reject SBSOFF-only runs.\n  const std::string globalcut_text = kin_info.globalcut.GetTitle();\n  if (globalcut_text.find("sbs.gem") != std::string::npos ||\n      globalcut_text.find("sbs.tr")  != std::string::npos) {\n    std::cerr << "ERROR: global_cut contains an SBS tracking requirement:\\n  "\n              << globalcut_text << "\\n"\n              << "Use the normal SBSOFF data config (for example GEN2_He3.cfg), "\n              << "not a config whose global_cut requires sbs.gem/sbs.tr." << std::endl;\n    return 3;\n  }\n\n  kin_info.rootfile_dir = sbsoff_dir;\n\n  std::cout << "SBSOFF primary data directory : " << sbsoff_dir << std::endl;\n  std::cout << "SBSON companion data directory: " << sbson_dir << std::endl;\n\n  getDB(kin_info.conf);\n\n  // Primary DATA chains.  is_data=1 is intentional: no simulation is loaded.\n  TChain *C = LoadRawRootFiles(kin_info, 1);\n  TChain *E = LoadRawRootFiles_E(kin_info, 1);\n\n  // SBSON is loaded one run at a time only when that run exists in the\n  // SBSOFF primary chain. This avoids duplicating physics events.\n  TChain *Csbs = nullptr;\n'''
    text = replace_once(text, old_load, new_load, "kinematic config / primary chain setup")

    # SBS variables are declared in the same place, but their input branches
    # must no longer be bound to the SBSOFF chain.
    text = replace_once(
        text,
        '  setrootvar::setbranch(C,"sbs.gem.track","nhits",&sbs_gem_nhits);\n',
        '  // sbs.gem.track is bound to the per-run SBSON companion chain below.\n',
        "SBS GEM branch binding",
    )
    text = replace_once(
        text,
        '  setrootvar::setbranch(C,"sbs.tr",sbstrvar,sbstrvar_mem);\n',
        '  // sbs.tr is bound to the per-run SBSON companion chain below.\n',
        "SBS track branch binding",
    )

    # Output event/matching flags.
    old_runbranch = '''  int T_runnum;         Tout->Branch("runnum", &T_runnum, "runnum/I");\n  TDatime T_datetime;   Tout->Branch("datetime", "TDatime", &T_datetime);\n'''
    new_runbranch = '''  int T_runnum;         Tout->Branch("runnum", &T_runnum, "runnum/I");\n  Long64_t T_evnum;     Tout->Branch("evnum", &T_evnum, "evnum/L");\n  bool T_run_has_sbson;     Tout->Branch("run_has_sbson", &T_run_has_sbson, "run_has_sbson/O");\n  bool T_sbson_event_match; Tout->Branch("sbson_event_match", &T_sbson_event_match, "sbson_event_match/O");\n  bool T_has_sbs_track;     Tout->Branch("has_sbs_track", &T_has_sbs_track, "has_sbs_track/O");\n  TDatime T_datetime;   Tout->Branch("datetime", "TDatime", &T_datetime);\n'''
    text = replace_once(text, old_runbranch, new_runbranch, "output run/event flags")

    # State for the currently loaded SBSON run.
    old_loopvars = '''  int treenum = 0, currenttreenum = 0, currentrunnum = 0;\n  int IHWP_run = -100;\n  time_t run_time_unix;\n'''
    new_loopvars = '''  int treenum = 0, currenttreenum = 0, currentrunnum = 0;\n  bool run_has_sbson = false;\n  int IHWP_run = -100;\n  time_t run_time_unix;\n'''
    text = replace_once(text, old_loopvars, new_loopvars, "event-loop state")

    # On a new run, construct a companion SBSON chain for only that run and
    # index it by CODA event number. Missing SBSON runs are allowed.
    #
    # This intentionally does NOT depend on the exact local spelling of the
    # filename/run-number parser. Locate the unique assignment to T_runnum
    # and insert the SBSON setup immediately after it.
    run_setup = r'''
      // Load/reload the SBSON companion only when the run number changes.
      if (T_runnum != currentrunnum) {
        currentrunnum = T_runnum;
        run_has_sbson = false;

        if (Csbs) {
          delete Csbs;
          Csbs = nullptr;
        }

        Utilities::KinConf sbs_run_info = kin_info;
        sbs_run_info.rootfile_dir = sbson_dir;
        sbs_run_info.runnums.clear();
        sbs_run_info.runnums.push_back(T_runnum);
        sbs_run_info.nruns = 1;

        // DATA only. If this run is absent from SBSON, the chain has 0 entries.
        Csbs = LoadRawRootFiles(sbs_run_info, 1);
        run_has_sbson = (Csbs && Csbs->GetEntries() > 0);

        if (run_has_sbson) {
          Csbs->SetBranchStatus("*", 0);
          Csbs->SetBranchStatus("g.evnum", 1);
          setrootvar::setbranch(Csbs,"sbs.gem.track","nhits",&sbs_gem_nhits);
          setrootvar::setbranch(Csbs,"sbs.tr",sbstrvar,sbstrvar_mem);

          const Long64_t nindexed = Csbs->BuildIndex("g.evnum");
          std::cout << "SBSON run " << T_runnum << ": "
                    << Csbs->GetEntries() << " entries, "
                    << nindexed << " indexed by g.evnum" << std::endl;
        } else {
          std::cout << "SBSON run " << T_runnum
                    << ": not available; keeping SBSOFF events with no SBS track."
                    << std::endl;
        }
      }

'''

    run_assign_re = re.compile(
        r'^(?P<indent>[ \t]*)T_runnum[ \t]*=[^;\n]+;[ \t]*$',
        re.MULTILINE,
    )
    run_matches = list(run_assign_re.finditer(text))

    if len(run_matches) != 1:
        nearby = "\n".join(
            line for line in text.splitlines()
            if "T_runnum" in line
        )
        raise RuntimeError(
            "Could not uniquely locate the T_runnum assignment for per-run "
            f"SBSON setup; found {len(run_matches)} candidates.\n"
            "Lines containing T_runnum were:\n" + nearby
        )

    m = run_matches[0]
    insert_at = m.end()
    text = text[:insert_at] + "\n" + run_setup + text[insert_at:]

    # Match the current SBSOFF event to the same CODA event in this run's SBSON replay.
    old_after_gcut = '''    if (!passedgCut) continue;\n\n    //cout<<"global passed "<<endl;\n'''
    new_after_gcut = '''    if (!passedgCut) continue;\n\n    // Match SBSOFF -> SBSON by CODA event number within the current run.\n    // Event numbers are not assumed to be globally unique across runs.\n    T_evnum = static_cast<Long64_t>(std::llround(evnum_T));\n    T_run_has_sbson = run_has_sbson;\n    T_sbson_event_match = false;\n    T_has_sbs_track = false;\n    ntrack_sbs = 0.0;\n\n    if (run_has_sbson && Csbs) {\n      const Long64_t sbson_entry = Csbs->GetEntryNumberWithIndex(T_evnum);\n      if (sbson_entry >= 0) {\n        Csbs->GetEntry(sbson_entry);\n        T_sbson_event_match = true;\n        T_has_sbs_track = (ntrack_sbs > 0.0);\n      }\n    }\n    \n    //cout<<"global passed "<<endl;\n'''
    text = replace_once(text, old_after_gcut, new_after_gcut, "event matching")

    # Guard the hadron angles against missing SBS tracks.
    old_angles = '''    double pphi = atan(py_sbs[0]/px_sbs[0]) ;\n    double ptheta = acos(pz_sbs[0]/p_sbs[0]);\n'''
    new_angles = '''    const double sbs_nan = std::numeric_limits<double>::quiet_NaN();\n    double pphi = sbs_nan;\n    double ptheta = sbs_nan;\n    if (T_has_sbs_track && p_sbs[0] > 0.0) {\n      pphi = std::atan2(py_sbs[0], px_sbs[0]);\n      const double cos_theta = std::max(-1.0, std::min(1.0, pz_sbs[0]/p_sbs[0]));\n      ptheta = std::acos(cos_theta);\n    }\n'''
    text = replace_once(text, old_angles, new_angles, "SBS angle guard")

    # Guard every output SBS-track quantity so SBSOFF-only runs / unmatched
    # events do not reuse stale array contents or read element [0] when absent.
    old_sbsout = '''    T_ntrack_sbs = ntrack_sbs;\n    T_ntrack_chi2_sbs = chi2_sbs[0];\n    T_ntrack_hits_sbs = sbs_gem_nhits[0];\n    T_vz_sbs = vz_sbs[0];\n    T_vx_sbs = vx_sbs[0];\n    T_vy_sbs = vy_sbs[0];\n    //T_xtgt_sbs = xtgt_sbs[0];\n    T_ytgt_sbs = ytgt_sbs[0];\n    T_thtgt_sbs = thtgt_sbs[0];\n    T_phtgt_sbs = phtgt_sbs[0];\n    //T_thetabend_sbs = thetabend_sbs;\n    T_xfp_sbs = xfp_sbs[0];\n    T_yfp_sbs = yfp_sbs[0];\n    T_thfp_sbs = thfp_sbs[0];\n    T_phfp_sbs = phfp_sbs[0];\n    T_trP_sbs = p_sbs[0];\n    T_trPx_sbs = px_sbs[0];\n    T_trPy_sbs = py_sbs[0];\n    T_trPz_sbs = pz_sbs[0];\n    T_trx_sbs = xTr_sbs[0];\n    T_try_sbs = yTr_sbs[0];\n    T_trth_sbs = thTr_sbs[0];\n    T_trph_sbs = phTr_sbs[0];\n'''
    new_sbsout = '''    T_ntrack_sbs = T_sbson_event_match ? ntrack_sbs : 0.0;\n    if (T_has_sbs_track) {\n      T_ntrack_chi2_sbs = chi2_sbs[0];\n      T_ntrack_hits_sbs = sbs_gem_nhits[0];\n      T_vz_sbs = vz_sbs[0];\n      T_vx_sbs = vx_sbs[0];\n      T_vy_sbs = vy_sbs[0];\n      //T_xtgt_sbs = xtgt_sbs[0];\n      T_ytgt_sbs = ytgt_sbs[0];\n      T_thtgt_sbs = thtgt_sbs[0];\n      T_phtgt_sbs = phtgt_sbs[0];\n      //T_thetabend_sbs = thetabend_sbs;\n      T_xfp_sbs = xfp_sbs[0];\n      T_yfp_sbs = yfp_sbs[0];\n      T_thfp_sbs = thfp_sbs[0];\n      T_phfp_sbs = phfp_sbs[0];\n      T_trP_sbs = p_sbs[0];\n      T_trPx_sbs = px_sbs[0];\n      T_trPy_sbs = py_sbs[0];\n      T_trPz_sbs = pz_sbs[0];\n      T_trx_sbs = xTr_sbs[0];\n      T_try_sbs = yTr_sbs[0];\n      T_trth_sbs = thTr_sbs[0];\n      T_trph_sbs = phTr_sbs[0];\n    } else {\n      T_ntrack_chi2_sbs = sbs_nan;\n      T_ntrack_hits_sbs = sbs_nan;\n      T_vz_sbs = sbs_nan;\n      T_vx_sbs = sbs_nan;\n      T_vy_sbs = sbs_nan;\n      T_ytgt_sbs = sbs_nan;\n      T_thtgt_sbs = sbs_nan;\n      T_phtgt_sbs = sbs_nan;\n      T_xfp_sbs = sbs_nan;\n      T_yfp_sbs = sbs_nan;\n      T_thfp_sbs = sbs_nan;\n      T_phfp_sbs = sbs_nan;\n      T_trP_sbs = sbs_nan;\n      T_trPx_sbs = sbs_nan;\n      T_trPy_sbs = sbs_nan;\n      T_trPz_sbs = sbs_nan;\n      T_trx_sbs = sbs_nan;\n      T_try_sbs = sbs_nan;\n      T_trth_sbs = sbs_nan;\n      T_trph_sbs = sbs_nan;\n    }\n'''
    text = replace_once(text, old_sbsout, new_sbsout, "guarded SBS output assignment")

    # Clean up the companion chain at the end.
    old_end = '''  fout->Write();\n  sw->Delete();\n\n  return 0;\n'''
    new_end = '''  fout->Write();\n  if (Csbs) {\n    delete Csbs;\n    Csbs = nullptr;\n  }\n  sw->Delete();\n  \n  return 0;\n'''
    text = replace_once(text, old_end, new_end, "SBSON cleanup")

    if text == original:
        raise RuntimeError("No changes were made")

    backup = macro.with_suffix(macro.suffix + ".before_dual_replay")
    if not backup.exists():
        shutil.copy2(macro, backup)

    macro.write_text(text, encoding="utf-8")

    print(f"Updated: {macro}")
    print(f"Backup : {backup}")
    print()
    print("Next checks:")
    print(f"  git -C {repo} diff -- scripts/replay/QuasiElastic_ana_withsbsgems.C")
    print("  # Then run ROOT with the normal SBSOFF config, not an SBS-track-cut config.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
