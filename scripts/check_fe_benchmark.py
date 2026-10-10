"""Compare the REAL Abaqus benchmark (runs/benchmark_griffith/raw) with the exact Griffith solution.

    python scripts/launch_abaqus.py --only benchmark_griffith --parallel 2
    python scripts/check_fe_benchmark.py

Checks, per crack length: G from Delta ALLSE vs exact incremental G (finite plate -> small
systematic offset of ~1-3 % is expected for W = 40a, H = 20a), crack-closure / energy ratio,
mode II share (should be ~0), and the FFM failure stress vs exact for several a/l_ch.
Pass criteria (suggested): |G error| < 3 % for da >= 2 elements, cc ratio within 0.97-1.03,
FFM strength error < 3 %. If this fails, fix the mesh/seam/NFORC extraction before the sweep.
"""
import sys, json, math, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.optimize import brentq
from ffm.kernel import Criterion
from ffm.materials import E_prime
from ffm.postprocess import Setup, evaluate

d = ROOT / "runs" / "benchmark_griffith"
st = Setup(d, raw_subdir="raw")
case = st.case
a = case["benchmark"]["a"]
mA = case["fe"]["materials"]["A"]
Ep = E_prime(mA["E"], mA["nu"], case["fe"]["plane"])
sig_ref = st.unc["F"] / case["W"]           # nominal stress at reference displacement
pd = st.path_data(case["paths"][0])
if pd is None:
    sys.exit("cracked benchmark results missing")

G_ex = math.pi * sig_ref ** 2 * (2 * a + pd.da) / (2 * Ep)
GI, GII = pd.modes()
print(f"E' = {Ep:.1f} MPa, sigma_ref = {sig_ref:.4f} MPa, {len(pd.da)} crack lengths")
if case["fe"].get("grips") == "sliding":
    s_uni = Ep * case["fe"]["delta"] / (2 * case["H"])
    print(f"sliding grips: sigma_ref / (E' * strain) = {sig_ref / s_uni:.4f}  (pre-crack compliance -> slightly < 1)")
else:
    print("WARNING: clamped grips - the centre stress is ~8 % above F/W; regenerate the benchmark case")
print("   da [mm]   G_dU/G_exact   closure/dU   G_II share")
for i in range(len(pd.da)):
    print(f"  {pd.da[i]:8.4f}   {pd.G_energy[i] / G_ex[i]:10.4f}   {pd.cc_ratio[i]:10.4f}   "
          f"{GII[i] / max(GI[i] + GII[i], 1e-300):8.4f}")

worst = 0.0
for sigma_c, G_c in ((50.0, 3.0), (50.0, 1.0), (50.0, 0.5)):   # a/l_ch ~ 0.24, 0.73, 1.5
    l_ch = Ep * G_c / sigma_c ** 2
    f = lambda b: 2 * l_ch / (math.pi * (a + b)) - (b * b - a * a) / (b * b)
    b = brentq(f, a * (1 + 1e-12), a + 1e3 * l_ch + 10 * a)
    s_ex = sigma_c * math.sqrt(1 - a * a / (b * b))
    r = evaluate(st, Criterion(sigma_c=sigma_c, G_Ic=G_c), {}, paths=[pd])[0]
    s_fe = r["lam_c"] * sig_ref
    flag = " (da-bound flag)" if r["at_lower"] or r["at_upper"] else ""
    if not flag:
        worst = max(worst, abs(s_fe / s_ex - 1))
    print(f"a/l_ch = {a / l_ch:6.3f}: FFM-FE {s_fe:7.3f} MPa, exact {s_ex:7.3f} MPa, "
          f"da_c {r['da_c']:.4f} vs {b - a:.4f}{flag}")
print(f"max strength error {worst:.2%}  ->", "PASSED" if worst < 0.03 else "CHECK MESH / EXTRACTION")
