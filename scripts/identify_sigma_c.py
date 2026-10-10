"""Inverse identification of the interface tensile strength sigma_c (PLAN.md, IP3).

A face-to-face tensile test does NOT measure sigma_c directly: the free-edge singularity at
the bimaterial interface lowers the failure stress. Given the measured nominal failure stress
and G_Ic, find the sigma_c for which the FFM prediction of the flat setup equals the test.

    python scripts/identify_sigma_c.py --pair PLA/PETG --sigma_f 18.5 --width 10
         [--G_Ic 0.6] [--sigma_f_sd 1.2]

--width : specimen width across the interface [mm]; the flat FE setup (width W) is scaled
          by kappa = width / W.  The flat setup has H/W = 2 with clamped grips; if the real
          specimen is very different, build a dedicated setup with make_cases (same workflow).
Prints sigma_c (and a +-1 SD band if --sigma_f_sd is given). tau_c is scaled with the
tau_c/sigma_c ratio in materials.json. If no solution exists (sigma_f above the
toughness-controlled limit), G_Ic is too low for the test -> re-check the G data.
"""
import sys, json, argparse, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scipy.optimize import brentq
from ffm.kernel import Criterion
from ffm.materials import load_materials
from ffm.postprocess import Setup, evaluate

ap = argparse.ArgumentParser()
ap.add_argument("--pair", required=True, help="e.g. PLA/PETG (as in materials.json)")
ap.add_argument("--sigma_f", type=float, required=True, help="measured nominal failure stress [MPa]")
ap.add_argument("--sigma_f_sd", type=float, default=0.0)
ap.add_argument("--width", type=float, required=True, help="specimen width [mm]")
ap.add_argument("--G_Ic", type=float, default=None)
a = ap.parse_args()

mats = load_materials(ROOT / "config" / "materials.json")
A, B = a.pair.split("/")
from ffm.materials import interface_props
I = interface_props(mats, A, B)
G_Ic = a.G_Ic if a.G_Ic is not None else I["G_Ic"]
ratio = (I["tau_c"] / I["sigma_c"]) if I.get("tau_c") else None
cands = sorted((ROOT / "runs").glob(f"flat__{A}-{B}__tension*"))
if not cands:
    sys.exit("run the flat setup first")
st = Setup(cands[0])
paths = [p for p in st.all_paths() if p.kind == "interface"]
kappa = a.width / st.case["W"]


def predicted(sc):
    c = Criterion(sigma_c=sc, tau_c=None if ratio is None else ratio * sc, G_Ic=G_Ic,
                  G_IIc=I.get("G_IIc"), energy_law=I.get("energy_law", "linear"))
    r = evaluate(st, c, {}, kappa=kappa, paths=paths)[0]
    return r["sigma_nom_c"], r


def solve(target):
    f = lambda sc: predicted(sc)[0] - target
    lo, hi = 0.5 * target, 20 * target
    if f(hi) < 0:
        return None
    return brentq(f, lo, hi, xtol=1e-4)


for tag, t in (("mean", a.sigma_f), ("-1SD", a.sigma_f - a.sigma_f_sd), ("+1SD", a.sigma_f + a.sigma_f_sd)):
    if tag != "mean" and a.sigma_f_sd == 0:
        continue
    sc = solve(t)
    if sc is None:
        print(f"{tag}: no solution - test strength exceeds the toughness-controlled limit for G_Ic = {G_Ic}")
        continue
    _, r = predicted(sc)
    flag = "  (da-bound flag: extend da range)" if r["at_lower"] or r["at_upper"] else ""
    print(f"{tag}: sigma_f = {t:.2f} MPa -> interface sigma_c = {sc:.2f} MPa "
          f"(onset {r['governing']}, da_c = {r['da_c']:.3f} mm, path {r['pid']}){flag}")
