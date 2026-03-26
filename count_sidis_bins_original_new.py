#!/usr/bin/env python3
"""
count_sidis_bins.py - The "Beast" Version
High-performance, memory-efficient SIDIS binning with Mean-Tracking.
"""

import argparse
import sys
import glob
import os
from itertools import product
from typing import Dict, List

import numpy as np
import pandas as pd

try:
    import boost_histogram as bh
except ImportError:
    sys.stderr.write("boost-histogram is required. Install with: python3 -m pip install --user boost-histogram\n")
    raise

# Extended default bins including nu
DEFAULT_BINS: Dict[str, List[float]] = {
    "Q2": [1.0, 2.0, 4.0, 8.0],
    "xB": [0.0, 0.2, 0.4, 0.6, 1.1],
    "nu": [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    "zh": [0.0, 0.3, 0.5, 0.7, 0.9, 1.1],
    "pT2": [0.0, 0.3, 0.7, 1.2, 2.0, 3.0],
    "phi_h": [0.0, 60.0, 120.0, 180.0, 240.0, 300.0, 360.0],
}

def parse_binspec(spec: str) -> Dict[str, List[float]]:
    out: Dict[str, List[float]] = {}
    if not spec: return out
    for part in [p.strip() for p in spec.split(";") if p.strip()]:
        if "=" not in part: continue
        dim, vals = part.split("=", 1)
        edges = [float(x) for x in vals.split(",") if x.strip()]
        if len(edges) >= 2: out[dim.strip()] = sorted(edges)
    return out

def ensure_bins(axis_list: List[str], user_bins: Dict[str, List[float]]) -> Dict[str, List[float]]:
    dims = set(axis_list)
    bins: Dict[str, List[float]] = {}
    for dim in dims:
        bins[dim] = user_bins.get(dim, DEFAULT_BINS.get(dim, [0, 1]))
    return bins

def build_axes(dim_order: List[str], bins: Dict[str, List[float]]):
    axes = []
    for dim in dim_order:
        edges = bins[dim]
        # Circular axis for phi_h if it spans 0-360
        if dim == "phi_h" and len(edges) >= 2 and edges[0] == 0.0 and edges[-1] == 360.0:
            axes.append(bh.axis.Regular(len(edges) - 1, edges[0], edges[-1], circular=True))
            continue
        axes.append(bh.axis.Variable(edges))
    return axes

def main():
    ap = argparse.ArgumentParser(description="Beast Mode SIDIS bin counter")
    ap.add_argument("parquet", nargs='+', help="Path to parquet files (wildcards supported)")
    ap.add_argument("--axes", default="zh", help="Comma list for hadron axes")
    ap.add_argument("--bins", help='Bin spec string e.g. "zh=0,0.2,0.4;pT2=0,1,2"')
    ap.add_argument("--apply-dis", action="store_true", help="Apply DIS cuts")
    ap.add_argument("--q2-min", type=float, default=1.0)
    ap.add_argument("--w-min", type=float, default=2.0)
    ap.add_argument("--y-min", type=float, default=0.25)
    ap.add_argument("--y-max", type=float, default=0.85)
    # Automated Hadron Cuts
    ap.add_argument("--zh-min", type=float, default=0.3)
    ap.add_argument("--zh-max", type=float, default=0.7)
    ap.add_argument("--pt2-max", type=float, default=1.2)

    ap.add_argument("--integrate-over", default="Q2,xB", help="Comma-separated dimensions to integrate")
    ap.add_argument("--out-csv", help="Output CSV path")
    args = ap.parse_args()

    file_list = []
    for p in args.parquet:
        file_list.extend(glob.glob(p))

    # Ensure list is unique and sorted for consistent RG-D analysis
    file_list = sorted(list(set(file_list)))    
    
    if not file_list:
        print(f"Error: No files found matching: {args.parquet}"); sys.exit(1)

    # Setup Axes
    axis_list = [a.strip() for a in args.axes.split(",") if a.strip()]
    integrate_list = [d.strip() for d in args.integrate_over.split(",") if d.strip()] if args.integrate_over else []
    integrate_set = set(integrate_list)
    
    dims_for_bins = list(set(axis_list) | integrate_set)
    bins = ensure_bins(dims_for_bins, parse_binspec(args.bins))
    
    # Initialize Histograms with Weight Storage
    hist_e = bh.Histogram(*build_axes(["Q2", "xB"], bins), storage=bh.storage.Weight())
    h_axes = build_axes(axis_list, bins)
    hist_h = bh.Histogram(*h_axes, storage=bh.storage.Weight())
    
    # Mean Trackers: Histograms to store sum(value) to calculate true mean later
    mean_trackers = {dim: bh.Histogram(*h_axes, storage=bh.storage.Double()) for dim in axis_list}

    # --- Helicity Yield Histograms ---
    hel_yields = {
        "N_plus":  bh.Histogram(*h_axes, storage=bh.storage.Weight()),
        "N_minus": bh.Histogram(*h_axes, storage=bh.storage.Weight()),
    }

    # pT2 moments (for pT-broadening): Σ w pT2 and Σ w pT2^2
    pt2_moments = {
        "pt2":    bh.Histogram(*h_axes, storage=bh.storage.Double()),
        "pt2_sq": bh.Histogram(*h_axes, storage=bh.storage.Double()),
    }
    
    print(f"Processing {len(file_list)} files iteratively...")
    for i, fpath in enumerate(file_list):
        df = pd.read_parquet(fpath)
        
        if args.apply_dis:
            e_mask = (df["w_e"] == 1) & (df["Q2"] >= args.q2_min) & (df["W"] >= args.w_min) & (df["y"] >= args.y_min) & (df["y"] <= args.y_max)
            passing_events = df.loc[e_mask, "sel_event_idx"].unique()
            df = df[df["sel_event_idx"].isin(passing_events)]

        df["hel_sign"] = np.where(
            df["helicity"] != 0,
            df["helicity"].astype(float),
            np.nan,
        )

        # Fill Electron side
        #df_e = df[df["w_e"] == 1].dropna(subset=["Q2", "xB"]) # OLD: blundered duplicate events    
        df_e = df[df["w_e"] == 1].drop_duplicates(subset=['run', 'sel_event_idx']).dropna(subset=["Q2", "xB"])
        hist_e.fill(df_e["Q2"].to_numpy(), df_e["xB"].to_numpy(), weight=df_e["w_e"].to_numpy())

        # Fill Hadron side
        df_h = df[
            (df["w_pip"] == 1)
            & (df["zh"] >= args.zh_min)
            & (df["zh"] <= args.zh_max)
            & (df["pT2"] <= args.pt2_max)
        ].dropna(subset=axis_list)

        if not df_h.empty:
            h_data = [df_h[dim].to_numpy() for dim in axis_list]
            h_weights = df_h["w_pip"].to_numpy()
            hist_h.fill(*h_data, weight=h_weights)
            
            # mean trackers for each axis dimension
            for dim in axis_list:
                mean_trackers[dim].fill(
                    *h_data,
                    weight=df_h[dim].to_numpy() * h_weights
                )

            # --- Fill Helicity Yields ---
            df_valid_hel = df_h.dropna(subset=["hel_sign"])
            if not df_valid_hel.empty:
                # Plus Helicity
                mask_p = df_valid_hel["hel_sign"] == 1.0
                if mask_p.any():
                    df_p = df_valid_hel[mask_p]
                    h_data_p = [df_p[dim].to_numpy() for dim in axis_list]
                    hel_yields["N_plus"].fill(*h_data_p, weight=df_p["w_pip"].to_numpy())
                
                # Minus Helicity
                mask_m = df_valid_hel["hel_sign"] == -1.0
                if mask_m.any():
                    df_m = df_valid_hel[mask_m]
                    h_data_m = [df_m[dim].to_numpy() for dim in axis_list]
                    hel_yields["N_minus"].fill(*h_data_m, weight=df_m["w_pip"].to_numpy())


            # pT2 moments for pT-broadening (even if pT2 is NOT in axes)
            if "pT2" in df_h.columns:
                pt2_vals = df_h["pT2"].to_numpy()
                pt2_moments["pt2"].fill(
                    *h_data,
                    weight=pt2_vals * h_weights
                )
                pt2_moments["pt2_sq"].fill(
                    *h_data,
                    weight=(pt2_vals ** 2) * h_weights
                )

        if (i+1) % 10 == 0: print(f"  - Progress: {i+1}/{len(file_list)} files")

# --- UNIVERSAL DYNAMIC NORMALIZATION (FIXED) ---
    view_e, var_e = hist_e.view().value, hist_e.view().variance
    view_h, var_h = hist_h.view().value, hist_h.view().variance
    axis_bins = [bins[d] for d in axis_list]
    rows = []
    
    # Map where Q2 (0) and xB (1) are in your hadron axes
    e_axes_in_hadron = {
        0: axis_list.index("Q2") if "Q2" in axis_list else None,
        1: axis_list.index("xB") if "xB" in axis_list else None
    }

    # Get views for helicity and pT2
    view_p = hel_yields["N_plus"].view()
    view_m = hel_yields["N_minus"].view()
    pt2_views = {name: hist.view() for name, hist in pt2_moments.items()}

    # --- INSERT THIS MISSING HELPER FUNCTION HERE ---
    def mean_and_err(sum_f, sum_f2, N):
        if N <= 0: return np.nan, np.nan
        mu = sum_f / N
        var = (sum_f2 / N) - mu * mu
        if var < 0: var = 0.0
        return mu, np.sqrt(var / N)
    # ------------------------------------------------

    print(f"\nGenerating final table with Universal Normalization...")
    for idx in product(*[range(len(b) - 1) for b in axis_bins]):
        
        # Start with full 2D electron counts
        cur_e, cur_v = view_e, var_e
        
        # STEP 1: Handle Q2 (axis 0 of hist_e)
        if e_axes_in_hadron[0] is not None:
            q_idx = idx[e_axes_in_hadron[0]]
            cur_e, cur_v = cur_e[q_idx, :], cur_v[q_idx, :]
        else:
            cur_e, cur_v = np.sum(cur_e, axis=0), np.sum(cur_v, axis=0)
            
        # STEP 2: Handle xB (axis 1 of hist_e)
        if e_axes_in_hadron[1] is not None:
            x_idx = idx[e_axes_in_hadron[1]]
            n_e, v_e = cur_e[x_idx], cur_v[x_idx]
        else:
            n_e, v_e = np.sum(cur_e), np.sum(cur_v)

        # STEP 3: Hadron Counts
        n_h, v_h = view_h[idx], var_h[idx]
        if n_e == 0: continue
        
        row = {}
        for i, dim in enumerate(axis_list):
            row[f"{dim}_lo"], row[f"{dim}_hi"] = axis_bins[i][idx[i]], axis_bins[i][idx[i]+1]
            row[f"{dim}_mean"] = mean_trackers[dim].view()[idx] / n_h if n_h > 0 else (row[f"{dim}_lo"] + row[f"{dim}_hi"])/2

        row.update({"N_e": n_e, "V_e": v_e, "N_pip": n_h, "V_pip": v_h })
    
        # --- Add Helicity Yields and Variances to Output ---
        row["N_plus"] = view_p[idx].value
        row["V_plus"] = view_p[idx].variance
        row["N_minus"] = view_m[idx].value
        row["V_minus"] = view_m[idx].variance
        # --- pT2 mean and error ---
        if n_h > 0:
            sum_pt2   = pt2_views["pt2"][idx]
            sum_pt2_2 = pt2_views["pt2_sq"][idx]
            mean_pt2, err_pt2 = mean_and_err(sum_pt2, sum_pt2_2, n_h)
            row["mean_pT2"] = mean_pt2
            row["err_pT2"]  = err_pt2
        else:
            row["mean_pT2"] = np.nan
            row["err_pT2"]  = np.nan
 
        rows.append(row)
    if args.out_csv:
        pd.DataFrame(rows).to_csv(args.out_csv, index=False)
        print(f"Successfully wrote {len(rows)} bins to {args.out_csv}")

if __name__ == "__main__":
    main()