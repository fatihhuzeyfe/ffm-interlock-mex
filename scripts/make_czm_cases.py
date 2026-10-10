"""Cohesive-zone (CZM) counterparts of selected FFM configurations (validation + post-onset).

    python scripts/make_czm_cases.py                 # default set (see SELECTION / SWEEP)
    python scripts/make_czm_cases.py --dry           # list only

For every selected configuration the FE case of the FFM study is re-used (same unit cell,
materials, raster, periodic boundary conditions, load) but
  * lengths are physical: the geometry (defined for w = 1) is scaled to the feature width w [mm]
    (a cohesive law has an intrinsic length, so the response is no longer size-independent);
  * the whole interface is a cohesive surface: quadratic traction initiation with
    k_sigma * (sigma_c, tau_c) and linear softening with k_G * G_Ic, G_IIc = ratio * G_Ic
    (power law with exponent 1 = the linear mixed-mode energy criterion used in the FFM);
  * no candidate crack paths (the crack path is free along the interface; bulk cracking is
    not part of the CZM model);
  * displacement control up to DELTA_FACTOR x the FFM onset displacement (until the whole
    interface has separated; afterwards only mechanical locking of the elastic bulk is left,
    whose limit - bulk fracture - is not part of the CZM model).
The interface is meshed with elements of at most 1/8 of the cohesive zone length estimate
lcz = 0.5 E* G_Ic / sigma_c^2 (and at most 0.05 w). Cases are written to runs_czm/<name>/case.json
and run with:  python scripts/launch_abaqus.py --runs runs_czm --parallel N
"""
import sys, json, argparse, pathlib
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ffm.materials import load_materials, interface_props

GEOMETRIES = ["flat", "sine_ar1", "U_ar1", "T_ar1", "dovetail_ar0.5_a10", "dovetail_ar1_a20", "dovetail_ar2_a30"]
PAIRS = ["PLA-PETG", "PLA-ABS"]
LOADS = ["tension", "shear"]
REF = dict(k_sigma=1.0, k_G=1.0, w=5.0, ratio=2.0, raster="r0")
SELECTION = [(f"{g}__{p}__{l}__{REF['raster']}", REF["k_sigma"], REF["k_G"], REF["w"], REF["ratio"])
             for g in GEOMETRIES for p in PAIRS for l in LOADS]
# master-curve check: lch/w varied through the interface toughness and the feature size
# (lch/w ~ 0.1 ... 6, strength- to toughness-controlled), interlocked and flat reference
SWEEP = [(f"{g}__PLA-PETG__tension__r0", 1.0, kg, w, 2.0)
         for g in ("dovetail_ar1_a20", "flat") for kg in (0.1, 0.25, 1.0) for w in (2.0, 5.0, 10.0)]
DELTA_FACTOR = 5.0          # far enough for the whole interface to separate


def scale_case(src, w):
    c = json.loads(json.dumps(src))
    sc = lambda xy: [xy[0] * w, xy[1] * w]
    c["W"], c["H"] = src["W"] * w, src["H"] * w
    c["interface_points"] = [sc(q) for q in src["interface_points"]]
    pts = np.asarray(c["interface_points"])
    S = np.r_[0.0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))]
    c["interface_breakpoints"] = [[float(x), float(y), float(s), "vertex"] for (x, y), s in zip(pts, S)]
    c["bulk_lines"], c["paths"] = [], []
    c.pop("precrack", None)
    c["geometry"] = dict(src["geometry"], w=w)
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--out", default="runs_czm")
    ap.add_argument("--results", default="results")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    mats = load_materials(ROOT / "config" / "materials.json")
    gov = pd.read_csv(ROOT / a.results / "governing.csv")
    out = ROOT / a.out
    n = 0
    seen = set()
    for setup, ks, kg, w, ratio in SELECTION + SWEEP:
        name = f"{setup}__ks{ks:g}__kg{kg:g}__w{w:g}"
        if name in seen:
            continue
        seen.add(name)
        src = json.loads((ROOT / a.runs / setup / "case.json").read_text())
        row = gov[(gov.setup == setup) & np.isclose(gov.k_sigma, ks) & np.isclose(gov.k_G, kg) &
                  np.isclose(gov.width_mm, w) & np.isclose(gov.GIIc_over_GIc, ratio)]
        if row.empty:
            print("   no FFM result for", name, "- skipped"); continue
        row = row.iloc[0]
        A, B = src["fe"]["materials"]["A"], src["fe"]["materials"]["B"]
        I = interface_props(mats, A["name"], B["name"])
        sig, tau, GI = I["sigma_c"] * ks, I["tau_c"] * ks, I["G_Ic"] * kg
        GII = ratio * GI
        Es = 2.0 / (1.0 / A["E"] + 1.0 / B["E"])
        lcz = 0.5 * Es * GI / sig ** 2
        h_if = min(0.05 * w, lcz / 8.0)
        c = scale_case(src, w)
        c["name"] = name
        c["project_root"] = str(ROOT)
        fe = c["fe"]
        fe["mesh"] = dict(fe["mesh"], h_global=fe["mesh"]["h_global"] * w, h_line=h_if, n_per_edge=2)
        delta_onset = float(row.lam_c) * src["fe"]["delta"] * w
        fe["delta"] = DELTA_FACTOR * delta_onset
        fe.pop("crack_contact", None)
        c["czm"] = dict(sigma_c=sig, tau_c=tau, G_Ic=GI, G_IIc=GII, mixed_mode="power", power=1.0,
                        K=10.0 * max(A["E"], B["E"]) / h_if, viscosity=1e-5,
                        inc0=0.001, inc_max=0.001, inc_min=1e-10, n_inc_max=50000,
                        lcz_estimate=lcz, h_interface=h_if, k_sigma=ks, k_G=kg, w=w, GIIc_over_GIc=ratio,
                        source_setup=setup,
                        ffm=dict(lam_c=float(row.lam_c), sigma_nom_c=float(row.sigma_nom_c),
                                 da_c=float(row.da_c), pid=row.pid, kind=row.kind,
                                 delta_onset=delta_onset, lch_over_w=float(1.0 / row.mu_interface)))
        n += 1
        print(f"{name}: h_if = {h_if:.3f} mm (lcz ~ {lcz:.2f} mm), FFM onset {row.sigma_nom_c:.2f} MPa")
        if not a.dry:
            d = out / name
            d.mkdir(parents=True, exist_ok=True)
            (d / "case.json").write_text(json.dumps(c))
    print(f"{n} CZM cases {'(dry run)' if a.dry else '-> ' + str(out)}")


if __name__ == "__main__":
    main()
