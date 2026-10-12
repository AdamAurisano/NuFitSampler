"""Convert the student's MCMC output (MCMCResults/vectorFile.root) to .npy so the
validation report can compare it with the new samplers. Needs PyROOT, so run it with
the ROOT-enabled python:   python3 python/convert_legacy.py [in.root] [out.npy]
Output columns follow the new convention: s12, s13, s23, dcp[-180,180), dm3l[1e-3], dm21[1e-5].
"""

import sys

import numpy as np
import ROOT

src = sys.argv[1] if len(sys.argv) > 1 else "MCMCResults/vectorFile.root"
dst = sys.argv[2] if len(sys.argv) > 2 else "results/legacy_samples.npy"
f = ROOT.TFile.Open(src)
t = f.Get("vecTree")
t.GetEntry(0)
v = {b: np.array(list(getattr(t, b))) for b in ("dm221", "dm231", "t12", "t13", "t23", "dcp")}
dm21, dm31 = v["dm221"], v["dm231"]
dm3l = np.where(dm31 > 0, dm31, dm31 - dm21)  # NuFIT: Dm31 for NO, Dm32 = Dm31 - Dm21 for IO
dcp = (v["dcp"] + 180.0) % 360.0 - 180.0
out = np.c_[v["t12"], v["t13"], v["t23"], dcp, dm3l * 1e3, dm21 * 1e5]
np.save(dst, out)
print(f"{len(out)} legacy samples -> {dst}")
