import numpy as np
from common_cuts import is_fd_pip

# ---------- 1) Vz Windows (Positive Pions) ----------
# Keep the name generic like in electron_cuts.py
VZ_WINDOWS_OB = {
    "LD2": (-20.0000, 5.0000),
    "CxC": (-10.5319, 5.0000),
    "Cu":  (-9.76024, -5.1676),
    "Sn":  (-4.6413,  5.0000),
}

VZ_WINDOWS_IB = {
    "LD2": (-15.0000, 5.0000),
    "CxC": (-10.28,   5.0000),
    "Cu":  (-10.69,  -6.50),
    "Sn":  (-6.13,    5.0000),
}

# ---------- 2) Helper Functions ----------

def vz_mask(vz, target, polarity="OB"):
    """Generic name matching electron_cuts.py"""
    windows = VZ_WINDOWS_OB if polarity == "OB" else VZ_WINDOWS_IB
    vz_min, vz_max = windows[target]
    return (vz >= vz_min) & (vz <= vz_max)

def chi2pid_mask(chi2pid, cut_val=5.0):
    """Refined Pion identification."""
    return np.abs(np.asarray(chi2pid)) < cut_val

def dvz_mask(pip_vz, ele_vz, cut_val=2.0):
    """Delta Vz matching."""
    return np.abs(np.asarray(pip_vz) - ele_vz) < cut_val

# ---------- 3) The Master Cutflow Function ----------

def pip_cutflow(df, ele_vz, target, polarity="OB"):
    """
    Master function name is specific ('pip_cutflow'), 
    but internal logic uses generic helper names.
    """
    pid     = df["pid"].to_numpy()
    status  = df["status"].to_numpy()
    chi2pid = df["chi2pid"].to_numpy()
    vz      = df["vz"].to_numpy()
    
    masks = {}

    # STEP 1: BASE (PID 211 + FD Status)
    masks["base"] = is_fd_pip(pid, status)
    N_base = int(np.sum(masks["base"]))

    # STEP 2: chi2pid
    masks["chi2pid"] = masks["base"] & chi2pid_mask(chi2pid)
    
    # STEP 3: Vz (Generic helper call)
    masks["vz"] = masks["chi2pid"] & vz_mask(vz, target, polarity)
    
    # STEP 4: dvz (Generic helper call)
    masks["dvz"] = masks["vz"] & dvz_mask(vz, ele_vz)
    
    masks["final"] = masks["dvz"]

    # Generate Stats (Consistent dictionary keys)
    order = ["base", "chi2pid", "vz", "final"]
    cutflow = {step: {"N": int(np.sum(masks[step])), 
                      "eff_base": 100.0 * np.sum(masks[step]) / N_base if N_base > 0 else 0.0} 
               for step in order}

    return masks["final"], cutflow, masks