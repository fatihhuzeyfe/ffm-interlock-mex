"""Export the data behind every paper figure as plain CSV files (figures/data/), so the figures
can be drawn and edited in MATLAB (figures/matlab/) or any other tool.

    python analysis/export_figure_data.py [--root .] [--out figures/data]

Inputs : results/governing.csv, results/ranking.csv, runs/*/case.json (+ raw results for the
         FFM example curves), docs/mesh_convergence/convergence_table.csv (optional)
Outputs: one CSV per figure panel; every file starts with a header row. Column meanings are
         documented in figures/README.md.
"""
import sys, json, argparse, pathlib, itertools
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

EXAMPLE_SETUP = "dovetail_ar1_a20__PLA-PETG__tension__r0"
EXAMPLE_K = (1.0, 0.1)        # (k_sigma, k_G): onset inside the path (coupled stress-energy)
EXAMPLE_W = 5.0
GEOMETRY_GALLERY = ["flat__PLA-PETG__tension__r0", "sine_ar1__PLA-PETG__tension__r0",
                    "U_ar1__PLA-PETG__tension__r0", "T_ar1__PLA-PETG__tension__r0",
                    "dovetail_ar1_a20__PLA-PETG__tension__r0"]
PATH_GALLERY = "dovetail_ar1_a20__PLA-PETG__shear__r0"     # no mirror pruning -> all paths


# ----------------------------------------------------------------------------- helpers
def point_at(pts, S):
    pts = np.asarray(pts, float)
    seg = np.hypot(*np.diff(pts, axis=0).T)
    cum = np.r_[0.0, np.cumsum(seg)]
    i = int(np.clip(np.searchsorted(cum, S) - 1, 0, len(seg) - 1))
    t = (S - cum[i]) / seg[i]
    return pts[i] + t * (pts[i + 1] - pts[i])


def sub_polyline(pts, S0, S1):
    pts = np.asarray(pts, float)
    seg = np.hypot(*np.diff(pts, axis=0).T)
    cum = np.r_[0.0, np.cumsum(seg)]
    lo, hi = min(S0, S1), max(S0, S1)
    inner = [p for p, c in zip(pts, cum) if lo < c < hi]
    return np.vstack([point_at(pts, lo)] + inner + [point_at(pts, hi)])


def mechanism_class(root, setup, pid, kind, material, cache={}):
    """interface base / interface on feature / interface wavy-flat / bulk A / bulk B."""
    if kind == "bulk":
        return "bulk_A" if material == "A" else "bulk_B"
    key = (setup, pid)
    if key not in cache:
        c = json.loads((root / "runs" / setup / "case.json").read_text())
        p = next(q for q in c["paths"] if q["pid"] == pid)
        S = p["S0"] if p["mode"] == "centred" else p["S0"] + p["direction"] * p["L"] / 2
        y = point_at(c["interface_points"], S)[1]
        shape = c["geometry"]["shape"]
        cache[key] = ("interface_wavy" if shape in ("flat", "sine") else
                      "interface_base" if abs(y) < 1e-6 else "interface_feature")
    return cache[key]


def load_tables(root):
    r = pd.read_csv(root / "results" / "ranking.csv")
    g = pd.read_csv(root / "results" / "governing.csv")
    k = ["setup", "k_sigma", "k_G", "width_mm", "GIIc_over_GIc"]
    d = r.merge(g[k + ["mu_interface", "pid", "kind", "material", "interpenetration", "at_upper"]],
                on=k, suffixes=("", "_g"))
    d["geometry"] = d.setup.str.split("__").str[0]
    d["raster"] = d.setup.str.split("__").str[-1]
    d["lch_over_w"] = 1.0 / d.mu_interface
    d["mechanism"] = [mechanism_class(root, s, p, kd, m)
                      for s, p, kd, m in zip(d.setup, d.pid, d.kind, d.material)]
    return d


# ----------------------------------------------------------------------------- figure 1
def fig1(root, out):
    rows = []
    for s in GEOMETRY_GALLERY:
        c = json.loads((root / "runs" / s / "case.json").read_text())
        for i, (x, y) in enumerate(c["interface_points"]):
            rows.append(dict(geometry=s.split("__")[0], point=i, x=x, y=y, W=c["W"], H=c["H"]))
    pd.DataFrame(rows).to_csv(out / "fig1a_interfaces.csv", index=False)

    c = json.loads((root / "runs" / PATH_GALLERY / "case.json").read_text())
    rows = []
    for p in c["paths"]:
        if p["kind"] == "interface":
            S1 = p["S0"] + (p["direction"] * p["L"] if p["mode"] == "edge" else p["L"])
            S0 = p["S0"] if p["mode"] == "edge" else p["S0"] - p["L"]
            xy = sub_polyline(c["interface_points"], S0, S1)
        else:
            o = np.asarray(p["origin"]); e = np.array([np.cos(np.radians(p["angle_deg"])),
                                                       np.sin(np.radians(p["angle_deg"]))])
            xy = np.vstack([o, o + p["L"] * e])
        for i, (x, y) in enumerate(xy):
            rows.append(dict(pid=p["pid"], kind=p["kind"], point=i, x=x, y=y))
    pd.DataFrame(rows).to_csv(out / "fig1b_candidate_paths.csv", index=False)
    pd.DataFrame([dict(point=i, x=x, y=y) for i, (x, y) in enumerate(c["interface_points"])]) \
        .to_csv(out / "fig1b_interface.csv", index=False)
    pd.DataFrame([dict(pid=p["pid"], kind=p["kind"], k=k, da=cr["da"]) for p in c["paths"]
                  for k, cr in enumerate(p["cracks"])]).to_csv(out / "fig1b_crack_lengths.csv", index=False)
    # every geometry of the study (interface polylines, unit-cell size)
    rows = []
    for cj in sorted((root / "runs").glob("*__PLA-PETG__tension__r0/case.json")):
        c2 = json.loads(cj.read_text())
        for i, (x, y) in enumerate(c2["interface_points"]):
            rows.append(dict(geometry=cj.parent.name.split("__")[0], shape=c2["geometry"]["shape"],
                             ar=c2["geometry"].get("ar"), alpha_deg=c2["geometry"].get("alpha_deg"),
                             point=i, x=x, y=y, W=c2["W"]))
    pd.DataFrame(rows).to_csv(out / "fig_all_geometries.csv", index=False)

    # FFM curves lambda_s(da), lambda_G(da) of the governing path of the example (k=1, w=5 mm)
    from ffm.postprocess import Setup
    from ffm.kernel import ffm_path
    from ffm.materials import load_materials
    from run_postprocess import criteria_for
    st = Setup(root / "runs" / EXAMPLE_SETUP)
    mats = load_materials(root / "config" / "materials.json")
    ci, cb, A, B = criteria_for(st.case, mats, 2.0)
    ci = ci.scaled(EXAMPLE_K[0], EXAMPLE_K[1])          # interface quality of the example
    w = EXAMPLE_W
    kappa = w / st.case["geometry"].get("w", 1.0)
    F = st.unc["F"] / st.case["W"]
    rows, summ = [], []
    for pd_ in st.all_paths():
        if pd_.kind != "interface":
            continue
        GI, GII = pd_.modes()
        br = [(r * kappa, sn, tau) for r, sn, tau in pd_.branches]
        res = ffm_path(br, pd_.da * kappa, GI * kappa, GII * kappa, ci, return_curves=True)
        summ.append(dict(pid=pd_.pid, sigma_c_onset=res.lam_c * F, da_c=res.da_c, governing=res.governing,
                         k_sigma=EXAMPLE_K[0], k_G=EXAMPLE_K[1], w_mm=w))
        for a, ls, lg in zip(res.curves["da"], res.curves["lam_s"], res.curves["lam_G"]):
            rows.append(dict(pid=pd_.pid, da_mm=a, sigma_stress=ls * F, sigma_energy=lg * F))
    pd.DataFrame(rows).to_csv(out / "fig1c_ffm_curves.csv", index=False)
    pd.DataFrame(summ).sort_values("sigma_c_onset").to_csv(out / "fig1c_ffm_onset.csv", index=False)


# ----------------------------------------------------------------------------- figure 2
def fig2(d, out):
    il = d[d["shape"] != "flat"]
    bins = np.logspace(np.log10(il.lch_over_w.min()), np.log10(il.lch_over_w.max()), 19)
    rows = []
    for (load, pair), s in il[il.GIIc_over_GIc == 2].groupby(["load", "pair"]):
        b = pd.cut(s.lch_over_w, bins)
        q = s.groupby(b, observed=True).rel_to_flat.quantile([0.1, 0.5, 0.9]).unstack()
        n = s.groupby(b, observed=True).size()
        for iv, row in q.iterrows():
            rows.append(dict(load=load, pair=pair, lch_over_w=np.sqrt(iv.left * iv.right),
                             gain_p10=row[0.1], gain_median=row[0.5], gain_p90=row[0.9], n=int(n[iv])))
    pd.DataFrame(rows).to_csv(out / "fig2a_master_curve_binned.csv", index=False)
    # every point (for scatter / density plots)
    il[["load", "pair", "geometry", "shape", "raster", "k_sigma", "k_G", "width_mm", "GIIc_over_GIc",
        "lch_over_w", "rel_to_flat", "sigma_nom_c", "mechanism"]] \
        .to_csv(out / "fig2b_all_points.csv", index=False)
    # family medians in the two asymptotic regimes
    il2 = il[il.GIIc_over_GIc == 2].copy()
    il2["regime"] = pd.cut(il2.lch_over_w, [0, 0.5, 2.0, 1e9], labels=["strength", "transition", "toughness"])
    il2.groupby(["load", "shape", "regime"], observed=True).rel_to_flat.median().rename("gain_median") \
        .reset_index().to_csv(out / "fig2c_family_regimes.csv", index=False)


# ----------------------------------------------------------------------------- figure 3
def fig3(d, out):
    il = d[~d["shape"].isin(["flat", "sine"])]
    s = il[(il.GIIc_over_GIc == 2) & (il.width_mm == 5)]
    t = s.groupby(["load", "k_sigma", "k_G"]).mechanism.value_counts(normalize=True).rename("share").reset_index()
    t.to_csv(out / "fig3a_mechanism_shares.csv", index=False)
    s = d[(d["shape"] != "flat") & (d.k_sigma == 1) & (d.k_G == 1) & (d.GIIc_over_GIc == 2)]
    s.groupby(["load", "pair", "shape", "width_mm"]).rel_to_flat.median().rename("gain_median") \
        .reset_index().to_csv(out / "fig3b_size_effect.csv", index=False)
    s = d[(d.k_sigma == 1) & (d.k_G == 1) & (d.GIIc_over_GIc == 2) & (d.width_mm == 5) & (d.raster == "r0")]
    s[["geometry", "shape", "pair", "load", "rel_to_flat"]].rename(columns={"rel_to_flat": "gain"}) \
        .sort_values(["geometry", "pair", "load"]).to_csv(out / "fig3c_gain_table.csv", index=False)
    k = ["geometry", "pair", "load", "k_sigma", "k_G", "width_mm", "GIIc_over_GIc"]
    p = d.pivot_table(index=k, columns="raster", values="sigma_nom_c")
    (p["r90"] / p["r0"]).rename("r90_over_r0").reset_index().to_csv(out / "fig3d_raster_ratio.csv", index=False)


# ----------------------------------------------------------------------------- figure 4
def sobol_exact(df, factors, y):
    V = df[y].var(ddof=0)
    out = []
    for f in factors:
        S1 = df.groupby(f)[y].mean().var(ddof=0) / V
        ST = 1 - df.groupby([x for x in factors if x != f])[y].mean().var(ddof=0) / V
        out.append(dict(factor=f, S1=S1, ST=ST))
    return pd.DataFrame(out)


def fig4(d, out):
    il = d[d["shape"] != "flat"].copy()
    il["ln_gain"] = np.log(il.rel_to_flat)
    F = ["geometry", "pair", "raster", "k_sigma", "k_G", "width_mm", "GIIc_over_GIc"]
    rows = []
    for load, s in il.groupby("load"):
        t = sobol_exact(s, F, "ln_gain"); t.insert(0, "load", load); rows.append(t)
    pd.concat(rows).to_csv(out / "fig4a_sobol.csv", index=False)

    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel as C
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import r2_score
    sub = il[(il["shape"] == "dovetail") & (il.pair == "PLA/PETG") & (il.raster == "r0") & (il.GIIc_over_GIc == 2)]
    maps, cv = [], []
    for load, s in sub.groupby("load"):
        X = np.c_[s.ar, s.alpha_deg / 10, np.log10(s.lch_over_w)]
        y = s.ln_gain.values
        kern = C(0.1) * RBF([0.5, 0.5, 0.3], (1e-2, 1e2)) + WhiteKernel(1e-3, (1e-6, 1e-1))
        pr = np.zeros_like(y)
        for tr, te in GroupKFold(6).split(X, y, s.geometry):
            pr[te] = GaussianProcessRegressor(kern, normalize_y=True, random_state=0).fit(X[tr], y[tr]).predict(X[te])
        cv.append(dict(load=load, model="GP dovetail PLA/PETG", R2_leave_one_geometry_out=r2_score(y, pr)))
        gp = GaussianProcessRegressor(kern, normalize_y=True, random_state=0).fit(X, y)
        A, L = np.meshgrid(np.linspace(0.5, 2, 61), np.linspace(10, 30, 61))
        for lw in (0.5, 3.0, 20.0):
            Z = np.exp(gp.predict(np.c_[A.ravel(), L.ravel() / 10, np.full(A.size, np.log10(lw))]))
            maps.append(pd.DataFrame(dict(load=load, lch_over_w=lw, ar=A.ravel(), alpha_deg=L.ravel(), gain=Z)))
    pd.concat(maps).to_csv(out / "fig4b_design_map.csv", index=False)

    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.model_selection import cross_val_predict
    X = pd.DataFrame(dict(log_lch_w=np.log10(il.lch_over_w), shear=(il.load == "shear").astype(int),
                          ar=il.ar.fillna(0), alpha=il.alpha_deg.fillna(0),
                          r90=(il.raster == "r90").astype(int), ratio=il.GIIc_over_GIc,
                          k_sigma=il.k_sigma, k_G=il.k_G, w=il.width_mm))
    for s_ in ("T", "U", "dovetail", "sine"):
        X["is_" + s_] = (il["shape"] == s_).astype(int)
    mats = {"PLA": 3071.0, "PETG": 1711.0, "ABS": 1699.0}
    X["E_mismatch"] = [(mats[a] - mats[b]) / (mats[a] + mats[b]) for a, b in (p.split("/") for p in il.pair)]
    m = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.06, max_leaf_nodes=31, random_state=0)
    pr = cross_val_predict(m, X, il.ln_gain, cv=GroupKFold(6), groups=il.geometry)
    cv.append(dict(load="both", model="gradient boosting, all inputs",
                   R2_leave_one_geometry_out=r2_score(il.ln_gain, pr)))
    pd.DataFrame(dict(gain_ffm=il.rel_to_flat.values, gain_surrogate=np.exp(pr), load=il.load.values)) \
        .to_csv(out / "fig4c_surrogate_parity.csv", index=False)
    pd.DataFrame(cv).to_csv(out / "fig4_cv_scores.csv", index=False)


def verification(root, out):
    """Griffith benchmark (parsed from docs/benchmark_check.txt) and mesh convergence table."""
    import re, shutil
    b = root / "docs" / "benchmark_check.txt"
    if b.exists():
        txt = b.read_text(errors="ignore")
        rows = [dict(da=float(a), G_ratio=float(g), closure_over_energy=float(c), GII_share=float(m))
                for a, g, c, m in re.findall(r"^\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s*$", txt, re.M)]
        pd.DataFrame(rows).to_csv(out / "figS1a_benchmark_G.csv", index=False)
        rows = [dict(a_over_lch=float(a), sigma_ffm_fe=float(f), sigma_exact=float(e), da_fe=float(d1), da_exact=float(d2))
                for a, f, e, d1, d2 in re.findall(r"a/l_ch =\s*([\d.]+): FFM-FE\s+([\d.]+) MPa, exact\s+([\d.]+) MPa, "
                                                  r"da_c ([\d.]+) vs ([\d.]+)", txt)]
        pd.DataFrame(rows).to_csv(out / "figS1b_benchmark_strength.csv", index=False)
    m = root / "docs" / "mesh_convergence" / "convergence_table.csv"
    if m.exists():
        shutil.copy2(m, out / "figS1c_mesh_convergence.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    root = pathlib.Path(a.root)
    out = pathlib.Path(a.out) if a.out else root / "figures" / "data"
    out.mkdir(parents=True, exist_ok=True)
    d = load_tables(root)
    fig1(root, out); print("figure 1 data")
    fig2(d, out); print("figure 2 data")
    fig3(d, out); print("figure 3 data")
    fig4(d, out); print("figure 4 data")
    verification(root, out); print("verification data")
    print("written to", out)


if __name__ == "__main__":
    main()
