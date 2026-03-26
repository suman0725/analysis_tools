#!/usr/bin/env python3
"""
run_single_sample.py - Non-Interactive Batch Version
Optimized for RG-D Targets: LD2, Carbon, Copper, and Tin.

MODES:
1. Full Mode (Default): Produces 'sidis_TARGET_RUN.parquet' for physics analysis.
   - Saves clean electrons + pions with e_ and pip_ prefixes.
   
2. Diagnostic Mode (--diag-only): Produces 'diag_e_...' and 'diag_pip_...'.
   - Saves ALL candidate particles with pass/fail flags and labeled prefixes.
"""
import sys
import os
import argparse
from pathlib import Path
import time
import numpy as np
import pandas as pd
import uproot

sys.path.append(os.path.join(os.path.dirname(__file__), 'scripts'))
# Physics & Bank Module Imports
from physics_constants import E_BEAM, M_ELECTRON
from physics import (
    get_Q2, get_W, get_y, get_xB, get_nu,
    get_zh, get_pt2, get_phih, get_theta, get_phi,
)
from electron_cuts import electron_cutflow
from pip_cuts import pip_cutflow
from bank_builders import REC_BRANCHES, MC_BRANCHES, build_per_particle_arrays
from pids import PID
from common_cuts import detect_polarity, is_fd, is_fd_ele, is_fd_pip

# =========================================================
# HELPER: PION DIAGNOSTIC LOGIC (Case 1 & 2)
# =========================================================
def save_pip_diagnostic(df_all, final_e_mask, target, base, out_dir):
    """
    Handles the extraction and saving of Case 1 vs Case 2 Pions.
    Merges electron vertex (e_vz) and momentum for both cases.
    """
    # 1. Filter for FD Pi+ Candidates
    pip_mask = is_fd_pip(df_all["pid"], df_all["status"])
    df_pip = df_all[pip_mask].copy()

    # 2. Get the 'Best Candidate' electron for every event
    e_cand_mask = is_fd_ele(df_all['pid'], df_all['status'])
    df_e_cands = df_all[e_cand_mask].copy()
    df_e_best = df_e_cands.sort_values(["event_idx_local", "p"], ascending=[True, False]).groupby("event_idx_local").head(1)
    
    # 3. MERGE: Attach electron vertex and momentum to every pion row
    df_pip = pd.merge(df_pip, df_e_best[["event_idx_local", "px", "py", "pz", "vz"]], 
                      on="event_idx_local", how="left", suffixes=("", "_e_ref"))

    # 4. Rename and organize columns
    df_pip = df_pip.rename(columns={
        "p": "pip_p", "vz": "pip_vz", "beta": "pip_beta", 
        "chi2pid": "pip_chi2pid", "status": "pip_status",
        "vz_e_ref": "e_vz", # Electron vertex stored for both cases
        "dc_x_r1": "pip_dc_x_r1", "dc_y_r1": "pip_dc_y_r1",
        "dc_x_r2": "pip_dc_x_r2", "dc_y_r2": "pip_dc_y_r2",
        "dc_x_r3": "pip_dc_x_r3", "dc_y_r3": "pip_dc_y_r3"
    })

    # 5. SIDIS Physics (using the reference electron components)
    df_pip["pip_zh"]   = get_zh(E_BEAM, df_pip["px_e_ref"], df_pip["py_e_ref"], df_pip["pz_e_ref"], df_pip["px"], df_pip["py"], df_pip["pz"])
    df_pip["pip_pt2"]  = get_pt2(E_BEAM, df_pip["px_e_ref"], df_pip["py_e_ref"], df_pip["pz_e_ref"], df_pip["px"], df_pip["py"], df_pip["pz"])
    df_pip["pip_phih"] = get_phih(E_BEAM, df_pip["px_e_ref"], df_pip["py_e_ref"], df_pip["pz_e_ref"], df_pip["px"], df_pip["py"], df_pip["pz"], True)
    
    # 6. Basic Kinematics
    df_pip["pip_theta"] = get_theta(df_pip["px"], df_pip["py"], df_pip["pz"], True)
    df_pip["pip_phi"]   = get_phi(df_pip["px"], df_pip["py"], True)

    # 7. Flag Case 2: event contains a Golden Electron
    clean_e_events = df_all[final_e_mask]["event_idx_local"].unique()
    df_pip["pass_e_cleanness"] = df_pip["event_idx_local"].isin(clean_e_events)

    # Clean up temporary momentum columns, keeping e_vz
    df_pip = df_pip.drop(columns=["px_e_ref", "py_e_ref", "pz_e_ref"])

    # 8. Save to Parquet
    final_path = out_dir / f"diag_pip_{target}_{base}.parquet"
    df_pip.to_parquet(final_path)
    print(f"[DIAG] Saved Pion Diagnostics (with e_vz): {final_path}")

# ==========================
# 1. Process one ROOT file
# ==========================

def process_file(path, target, out_dir_path, sample_type="data",
                 max_events=None, apply_tm=False, tm_min=None,
                 start_event_idx=0, diag_only=False):
    
    pol = detect_polarity(path)
    branch_list = list(REC_BRANCHES)
    if sample_type == "sim" and apply_tm:
        branch_list += MC_BRANCHES

    with uproot.open(path) as f:
        # Robust Tree Name Check
        tree_name = "data" 
        if "data" not in f: 
            if "clas12" in f: tree_name = "clas12"
            elif "events" in f: tree_name = "events"
            
        tree = f[tree_name]
        arrs = tree.arrays(branch_list, library="ak", entry_stop=max_events)

    df_all = build_per_particle_arrays(arrs, target_group=target)

    if sample_type == "sim" and apply_tm:
        from truth_matching import add_truth_matching
        df_all = add_truth_matching(df_all, arrs, quality_min=tm_min)

    # --- THE ACCOUNTANT (Inventory Report) ---
    f_mask = is_fd(df_all["status"].to_numpy())
    df_fd = df_all[f_mask]
    print(f"\n{'='*60}\n ACCOUNTANT REPORT: {os.path.basename(path)}\n{'-'*60}")
    
    c_pos, c_neg = df_all['charge'] > 0, df_all['charge'] < 0
    fd_pos, fd_neg = df_fd['charge'] > 0, df_fd['charge'] < 0

    counts = {
        "Total Tracks":   (len(df_all), len(df_fd)),
        "Positives (+)":  (len(df_all[c_pos]), len(df_fd[fd_pos])),
        "Negatives (-)":  (len(df_all[c_neg]), len(df_fd[fd_neg])),
        "Electrons (e-)": (len(df_all[df_all['pid'] == PID.ELECTRON]),  len(df_fd[df_fd['pid'] == PID.ELECTRON])),
        "Pions (pi+)":    (len(df_all[df_all['pid'] == PID.PION_PLUS]), len(df_fd[df_fd['pid'] == PID.PION_PLUS])),
    }
    for lbl, (ac, fc) in counts.items():
        print(f" {lbl:<19} | {ac:<13} | {fc:<12}")
    print(f"{'='*60}\n")

    # --- ELECTRON SELECTION ---
    final_e_mask, e_cf, e_masks = electron_cutflow(df_all, target, pol, sample_type)

    # =========================================================
    # MODE 1: DIAGNOSTIC SAVING
    # =========================================================
    if diag_only:
        base = os.path.splitext(os.path.basename(path))[0]
        
        # 1. Electron Diagnostic with e_ Prefix
        e_cand_mask = is_fd_ele(df_all['pid'], df_all['status'])
        df_e_diag = df_all[e_cand_mask].copy()
        
        df_e_diag = df_e_diag.rename(columns={
            "p": "e_p", "vz": "e_vz", "sf": "e_sf", "phi": "e_phi",
            "dc_x_r1": "e_dc_x_r1", "dc_y_r1": "e_dc_y_r1",
            "dc_x_r2": "e_dc_x_r2", "dc_y_r2": "e_dc_y_r2",
            "dc_x_r3": "e_dc_x_r3", "dc_y_r3": "e_dc_y_r3"
        })
        
        px, py, pz = df_e_diag["px"].to_numpy(), df_e_diag["py"].to_numpy(), df_e_diag["pz"].to_numpy()
        df_e_diag["e_Q2"] = get_Q2(E_BEAM, px, py, pz)
        df_e_diag["e_xB"] = get_xB(E_BEAM, px, py, pz)
        df_e_diag["e_W"]  = get_W(E_BEAM, px, py, pz)
        df_e_diag["e_y"]  = get_y(E_BEAM, px, py, pz)
        df_e_diag["e_nu"] = get_nu(E_BEAM, px, py, pz)
        df_e_diag["e_theta"] = get_theta(px, py, pz, True)
        
        for step, m_arr in e_masks.items():
            df_e_diag[f"pass_{step}"] = m_arr[e_cand_mask]
            
        df_e_diag.to_parquet(out_dir_path / f"diag_e_{target}_{base}.parquet")
        save_pip_diagnostic(df_all, final_e_mask, target, base, out_dir_path)
        
        return pd.DataFrame(), e_cf, {}, len(arrs), 0

    # =========================================================
    # MODE 2: FINALIZED FULL SIDIS PRODUCTION
    # =========================================================
    df_e_all = df_all[final_e_mask].copy()
    if df_e_all.empty: return pd.DataFrame(), e_cf, {}, len(arrs), 0

    # Pick the best electron for DIS normalization
    df_e_top = df_e_all.sort_values(["event_idx_local", "p"], ascending=[True, False]).groupby("event_idx_local", as_index=False).head(1)
    
    # --- PION SELECTION & STATS ---
    # Only consider pions in events that passed the electron cuts
    df_pip_cands = df_all[df_all["event_idx_local"].isin(df_e_top["event_idx_local"].unique())].copy()
    
    # Merge electron vz into pion candidates for Delta-Vz matching
    df_pip_cands = pd.merge(df_pip_cands, df_e_top[["event_idx_local", "vz"]], 
                            on="event_idx_local", how="left", suffixes=("", "_e"))

    final_pip_mask, pip_cf, pip_masks = pip_cutflow(df_pip_cands, df_pip_cands["vz_e"], target, pol)

    # Isolate strictly final pions
    df_pip_final = df_pip_cands[final_pip_mask].copy()
    
    # Rename for clear e_ vs pip_ prefixes
    df_e_m = df_e_top.rename(columns={"px":"e_px", "py":"e_py", "pz":"e_pz", "p":"e_p", "vz":"e_vz", "event_id":"rc_event"})
    
    # Calculate DIS Kinematics
    px_e, py_e, pz_e = df_e_m["e_px"].to_numpy(), df_e_m["e_py"].to_numpy(), df_e_m["e_pz"].to_numpy()
    df_e_m["Q2"] = get_Q2(E_BEAM, px_e, py_e, pz_e)
    df_e_m["xB"] = get_xB(E_BEAM, px_e, py_e, pz_e)
    df_e_m["W"]  = get_W(E_BEAM, px_e, py_e, pz_e)
    df_e_m["nu"] = get_nu(E_BEAM, px_e, py_e, pz_e)
    df_e_m["w_e"] = 1  # Standard weight for DIS normalization
    df_e_m["y"] = get_y(E_BEAM, px_e, py_e, pz_e)

    df_pip_m = df_pip_final.rename(columns={"px":"pip_px", "py":"pip_py", "pz":"pip_pz", "p":"pip_p", "vz":"pip_vz", "beta":"pip_beta", "chi2pid":"pip_chi2pid"})

    # LEFT MERGE: Preserve all electrons (even if 0 pions) for RA denominator
    df_sidis = pd.merge(df_e_m, df_pip_m[["event_idx_local", "pip_px", "pip_py", "pip_pz", "pip_p", "pip_vz", "pip_beta", "pip_chi2pid"]], on="event_idx_local", how="left")
    
    if not df_sidis.empty:
        # SIDIS Variables (pT2 and zh are used for RA and Broadening)
        df_sidis["zh"] = get_zh(E_BEAM, df_sidis["e_px"], df_sidis["e_py"], df_sidis["e_pz"], df_sidis["pip_px"], df_sidis["pip_py"], df_sidis["pip_pz"])
        df_sidis["pT2"] = get_pt2(E_BEAM, df_sidis["e_px"], df_sidis["e_py"], df_sidis["e_pz"], df_sidis["pip_px"], df_sidis["pip_py"], df_sidis["pip_pz"])
        # phi_h is used for Azimuthal Modulations
        df_sidis["phi_h"] = get_phih(E_BEAM, df_sidis["e_px"], df_sidis["e_py"], df_sidis["e_pz"], df_sidis["pip_px"], df_sidis["pip_py"], df_sidis["pip_pz"], True)
        
        # Flag rows with a matched pion
        df_sidis["w_pip"] = np.where(df_sidis["pip_p"].notna(), 1, 0)
        df_sidis["sel_event_idx"] = df_sidis["event_idx_local"] + start_event_idx

    # Final columns needed for RA, Broadening, and Modulations
    out_cols = ["run", "rc_event", "sel_event_idx", "helicity", "helicityRaw", "w_e", "w_pip", "e_p", "e_vz", "Q2", "xB", "W", "nu", "y", 
                "pip_p", "pip_vz", "pip_beta", "pip_chi2pid", "zh", "pT2", "phi_h"]
    
    for c in out_cols:
        if c not in df_sidis.columns: df_sidis[c] = np.nan
            
    return df_sidis[out_cols], e_cf, pip_cf, len(arrs), len(df_e_top)

# ==========================
# 2. Process Wrapper
# ==========================

def process_target(files, target, out_dir_path, sample_type="data", max_events=None, 
                   diag_only=False, apply_tm=False, tm_min=0.4): 
    all_rows = []
    e_total_cf = None
    pip_total_cf = None
    next_idx = 0
    if not files: return pd.DataFrame(), None, None, "OB"
    pol = detect_polarity(files[0])

    for path in files:
        print(f"  - processing {path}")
        try:
            df, ecf, pcf, n_read, n_e = process_file(path, target, out_dir_path, sample_type, max_events, apply_tm, tm_min, next_idx, diag_only)
            
            # Electron cutflow aggregate
            if e_total_cf is None: e_total_cf = {k: {"N": 0} for k in ecf.keys()}
            for k in ecf: e_total_cf[k]["N"] += ecf[k]["N"]
            
            # Pion cutflow aggregate
            if pcf:
                if pip_total_cf is None: pip_total_cf = {k: {"N": 0} for k in pcf.keys()}
                for k in pcf: pip_total_cf[k]["N"] += pcf[k]["N"]

            if diag_only: continue
            next_idx += n_e
            if not df.empty: all_rows.append(df)
        except Exception as e:
            print(f"    !!! Error processing {path}: {e}")
            continue

    # Efficiencies
    if e_total_cf:
        base_n = e_total_cf.get("base", {}).get("N", 0)
        for k in e_total_cf: e_total_cf[k]["eff_base"] = 100.0 * e_total_cf[k]["N"] / base_n if base_n > 0 else 0.0
    
    if pip_total_cf:
        base_n = pip_total_cf.get("base", {}).get("N", 0)
        for k in pip_total_cf: pip_total_cf[k]["eff_base"] = 100.0 * pip_total_cf[k]["N"] / base_n if base_n > 0 else 0.0

    return pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame(), e_total_cf, pip_total_cf, pol

# ==========================
# 3. Batch Main Execution
# ==========================

def main():
    t_start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True); parser.add_argument("--root-file", required=True)
    parser.add_argument("--sim", action="store_true"); parser.add_argument("--entry-stop", type=int)
    parser.add_argument("--out-dir", default="."); parser.add_argument("--diag-only", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    df_sidis, e_cf, pip_cf, pol = process_target([args.root_file], args.target, out_dir, "sim" if args.sim else "data", args.entry_stop, args.diag_only)

    if not args.diag_only and not df_sidis.empty:
        out_name = f"sidis_{args.target}_{pol}_{os.path.basename(args.root_file)}.parquet"
        df_sidis.to_parquet(out_dir / out_name)

    # PRINT ELECTRON REPORT
    if e_cf:
        print(f"\n[ELECTRON CUTFLOW: {args.target}]")
        for k in ["base", "p", "sf", "pcal", "dc", "vz", "final"]:
            if k in e_cf: print(f"  {k:<12} | {e_cf[k]['N']:<10} | {e_cf[k]['eff_base']:>10.2f}%")

    # PRINT PION REPORT
    if pip_cf:
        print(f"\n[PION CUTFLOW: {args.target}]")
        for k in ["base", "chi2pid", "vz", "dvz", "final"]:
            if k in pip_cf: print(f"  {k:<12} | {pip_cf[k]['N']:<10} | {pip_cf[k]['eff_base']:>10.2f}%")

    print(f"\nDone in {time.time()-t_start:.1f}s")

if __name__ == "__main__": main()