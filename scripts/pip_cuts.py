import numpy as np
from common_cuts import is_fd_pip

# ---------- 1) Vz Windows (Positive Pions) ----------
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

# ---------- 2) Delta Vz Windows (pip_vz - ele_vz) ----------
# just loose cut, need to be updated 
DVZ_WINDOWS_OB = {
    "LD2": (-20.000, 20.0000),  
    "CxC": (-20.000, 20.0000), 
    "Cu":  (-20.000, 20.0000),
    "Sn":  (-20.000, 20.0000),
}

DVZ_WINDOWS_IB = {
    "LD2": (-20.000, 20.0000),
    "CxC": (-20.000, 20.0000),
    "Cu":  (-20.000, 20.0000),
    "Sn":  (-20.000, 20.0000),
}

# ---------- 3) Helper Functions ----------

def vz_mask(vz, target, polarity="OB"):
    """Vertex Z cut based on target and polarity."""
    windows = VZ_WINDOWS_OB if polarity == "OB" else VZ_WINDOWS_IB
    vz_min, vz_max = windows[target]
    return (vz >= vz_min) & (vz <= vz_max)

def dvz_mask(pip_vz, ele_vz, target, polarity="OB"):
    """Delta Vz matching using target-specific dictionaries."""
    windows = DVZ_WINDOWS_OB if polarity == "OB" else DVZ_WINDOWS_IB
    dvz_min, dvz_max = windows[target]
    dvz = np.asarray(pip_vz) - ele_vz
    return (dvz >= dvz_min) & (dvz <= dvz_max)

def chi2pid_mask(chi2pid, cut_val=10.0):
    """Refined Pion identification."""
    return np.abs(np.asarray(chi2pid)) < cut_val

# ---------- 4) The Master Cutflow Function ----------

def pip_cutflow(df, ele_vz, target, polarity="OB"):
    """
    Master function to identify pions and calculate cutflow stats.
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
    
    # STEP 4: dvz (Now uses the dictionary values from the top)
    masks["dvz"] = masks["vz"] & dvz_mask(vz, ele_vz, target, polarity)
    
    masks["final"] = masks["dvz"]

    # Generate Stats (Consistent dictionary keys)
    order = ["base", "chi2pid", "vz", "final"]
    cutflow = {step: {"N": int(np.sum(masks[step])), 
                      "eff_base": 100.0 * np.sum(masks[step]) / N_base if N_base > 0 else 0.0} 
               for step in order}

    return masks["final"], cutflow, masks