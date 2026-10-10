"""FFM post-processing of all finished FE setups.

    python scripts/run_postprocess.py                       # all setups in runs/ with raw/uncracked.json
    python scripts/run_postprocess.py --only dovetail       # subset (substring match)
    python scripts/run_postprocess.py --plots               # + FFM curves of every setup (k=1, middle size)
    python scripts/run_postprocess.py --benchmark           # synthetic Griffith check of this chain
    python scripts/run_postprocess.py --jobs 8              # parallel over setups

Strengths/toughnesses, interface quality (k_sigma, k_G) and physical feature width
enter ONLY here, so everything can be re-evaluated in minutes without re-running Abaqus.

Outputs (results/):
    paths.csv       every path x k_sigma x k_G x width x mode-ratio   (onset load of each candidate path)
    governing.csv   the weakest path per setup/parameter combination   (= predicted onset location and load)
    ranking.csv     governing sigma_nom_c pivoted by geometry, normalised by the flat interface
    curves/*.png    lambda_s, lambda_G vs crack length (with --plots)
"""
import sys, json, argparse, pathlib, itertools, csv, traceback
from concurrent.futures import ProcessPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from ffm.kernel import Criterion, ffm_path
from ffm.materials import load_materials, E_star, E_prime, brittleness, bulk_props_for_crack, interface_props
from ffm.postprocess import Setup, evaluate


# ----------------------------------------------------------------------------- criteria
def criteria_for(case, mats, gii_ratio=None):
    """Interface criterion and bulk criteria of one setup (reference quality k = 1).
    Orthotropic materials: the bulk criterion depends on the crack direction relative to the
    roads -> returned as a function angle_deg -> Criterion (see ffm.materials)."""
    fA, fB = case["fe"]["materials"]["A"], case["fe"]["materials"]["B"]
    A, B = fA["name"], fB["name"]
    M = mats["materials"]

    def bulk(fe_mat):
        m = M[fe_mat["name"]]
        raster = fe_mat.get("raster_deg", 0.0)
        if "X" in m and fe_mat.get("ortho"):
            def f(angle_deg, m=m, raster=raster):
                sc, gc = bulk_props_for_crack(m, angle_deg, raster)
                return Criterion(sigma_c=sc, G_Ic=gc, energy_law="total")
            return f
        return Criterion(sigma_c=m.get("sigma_c", m.get("Y")), G_Ic=m.get("G_c", m.get("G_L")),
                         energy_law="total")
    cb = {"A": bulk(fA), "B": bulk(fB)}
    if A == B and f"{A}/{B}" not in mats["interfaces"]:
        b = cb["A"]
        return (b(0.0) if callable(b) else b), cb, A, B     # homogeneous (benchmark)
    I = interface_props(mats, A, B)
    G_IIc = I.get("G_IIc")
    if gii_ratio is not None:
        G_IIc = gii_ratio * I["G_Ic"]
    ci = Criterion(sigma_c=I["sigma_c"], tau_c=I.get("tau_c"), G_Ic=I["G_Ic"], G_IIc=G_IIc,
                   energy_law=I.get("energy_law", "linear"), power=I.get("power", 1.0),
                   bk_eta=I.get("bk_eta", 2.0))
    return ci, cb, A, B


# ----------------------------------------------------------------------------- one setup
def process_setup(args):
    d, study, mats, raw_subdir, plots, plot_dir = args
    try:
        st = Setup(d, raw_subdir=raw_subdir)
    except FileNotFoundError:
        return d.name, [], "no uncracked result"
    paths = st.all_paths()
    if not paths:
        return d.name, [], "no cracked results"
    case = st.case
    n_expected = len(case["paths"])
    post = study["post"]
    ks, kg = post.get("k_sigma", [1.0]), post.get("k_G", [1.0])
    widths = post.get("feature_width_mm", [case["geometry"].get("w", 1.0)])
    ratios = post.get("GIIc_over_GIc", [None])
    w_fe = case["geometry"].get("w", 1.0)
    plane = case["fe"]["plane"]
    mA, mB = case["fe"]["materials"]["A"], case["fe"]["materials"]["B"]
    Es = E_star(mA["E"], mA["nu"], mB["E"], mB["nu"], plane) if mA["name"] != mB["name"] \
        else E_prime(mA["E"], mA["nu"], plane)
    g = case["geometry"]
    rows = []
    for ratio in ratios:
        ci0, cb, A, B = criteria_for(case, mats, ratio)
        for k_s, k_g, w in itertools.product(ks, kg, widths):
            ci = ci0.scaled(k_s, k_g)
            kappa = w / w_fe
            res = evaluate(st, ci, cb, kappa=kappa, paths=paths)
            mu = brittleness(w, Es, ci.G_Ic, ci.sigma_c)          # sigma_c^2 w / (E* G_Ic)
            for r in res:
                r.update(shape=g["shape"], ar=g.get("ar"), alpha_deg=g.get("alpha_deg"),
                         pair=f"{A}/{B}", load=case["fe"]["load"], k_sigma=k_s, k_G=k_g,
                         width_mm=w, GIIc_over_GIc=ratio, mu_interface=mu,
                         paths_done=len(paths), paths_expected=n_expected)
            rows.extend(res)
    if plots:
        _plot_curves(st, paths, case, mats, widths, w_fe, plot_dir)
    return d.name, rows, ""


def _plot_curves(st, paths, case, mats, widths, w_fe, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ci, cb, A, B = criteria_for(case, mats)
    w = widths[len(widths) // 2]
    kappa = w / w_fe
    fig, ax = plt.subplots(figsize=(7, 4.5))
    best = None
    for pd in paths:
        crit = ci if pd.kind == "interface" else cb[pd.material]
        if callable(crit):
            crit = crit(pd.meta.get("angle_deg") or 0.0)
        GI, GII = pd.modes()
        br = [(r * kappa, sn, tau) for r, sn, tau in pd.branches]
        res = ffm_path(br, pd.da * kappa, GI * kappa, GII * kappa, crit, return_curves=True)
        c = res.curves
        l, = ax.loglog(c["da"], c["lam_s"], "-", lw=1, alpha=.8)
        ax.loglog(c["da"], c["lam_G"], "--", lw=1, color=l.get_color(), alpha=.8)
        ax.plot(res.da_c, res.lam_c, "o", color=l.get_color(), ms=4)
        if best is None or res.lam_c < best[0]:
            best = (res.lam_c, res.da_c, pd.pid)
    ax.plot(best[1], best[0], "k*", ms=12, label=f"governing: {best[2]}")
    ax.set_ylim(0.5 * best[0], 30 * best[0])     # paths blocked by compression -> lambda = inf
    ax.set_xlabel("crack increment $\\Delta a$ [mm]")
    ax.set_ylabel("load factor $\\lambda$  (solid: stress, dashed: energy)")
    ax.set_title(f"{case['name']}  (w = {w} mm, k = 1)", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / f"{case['name']}.png", dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------- tables
KEYS_PATH = ["setup", "shape", "ar", "alpha_deg", "pair", "load", "k_sigma", "k_G", "width_mm",
             "GIIc_over_GIc", "mu_interface", "pid", "kind", "material", "note", "lam_c",
             "sigma_nom_c", "da_c", "governing", "at_lower", "at_upper", "segment_exhausted", "cc_ratio_median",
             "interpenetration", "overlap_ratio", "pen_share", "lam_c_nopen", "n_closure_only", "mode_mix_at_c", "paths_done", "paths_expected"]
GROUP = ["setup", "k_sigma", "k_G", "width_mm", "GIIc_over_GIc"]


def write_csv(path, rows, keys):
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        for r in rows:
            wr.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in keys})


def governing_rows(rows):
    best = {}
    second = {}
    for r in rows:
        key = tuple(r[k] for k in GROUP)
        if key not in best or r["lam_c"] < best[key]["lam_c"]:
            if key in best:
                second[key] = best[key]
            best[key] = r
        elif key not in second or r["lam_c"] < second[key]["lam_c"]:
            second[key] = r
    out = []
    for key, r in best.items():
        g = dict(r)
        s = second.get(key)
        # competition margin: how close the runner-up is (small -> location sensitive to data)
        g["runner_up"] = s["pid"] if s else ""
        g["runner_up_kind"] = s["kind"] if s else ""
        g["margin"] = (s["lam_c"] / r["lam_c"] - 1) if s else float("nan")
        g["complete"] = r["paths_done"] == r["paths_expected"]
        out.append(g)
    return out


def ranking_rows(gov):
    """sigma_nom_c of each geometry normalised by the flat interface of the same pair/load/params."""
    flat = {}
    for g in gov:
        if g["shape"] == "flat":
            flat[(g["pair"], g["load"], g["k_sigma"], g["k_G"], g["width_mm"], g["GIIc_over_GIc"])] = g["sigma_nom_c"]
    out = []
    for g in gov:
        key = (g["pair"], g["load"], g["k_sigma"], g["k_G"], g["width_mm"], g["GIIc_over_GIc"])
        f = flat.get(key)
        out.append(dict(pair=g["pair"], load=g["load"], k_sigma=g["k_sigma"], k_G=g["k_G"],
                        width_mm=g["width_mm"], GIIc_over_GIc=g["GIIc_over_GIc"], setup=g["setup"],
                        shape=g["shape"], ar=g["ar"], alpha_deg=g["alpha_deg"],
                        sigma_nom_c=g["sigma_nom_c"], rel_to_flat=(g["sigma_nom_c"] / f) if f else "",
                        onset_location=f"{g['kind']}:{g['pid']} ({g['note']})",
                        governing=g["governing"], margin=g["margin"]))
    out.sort(key=lambda r: (r["pair"], r["load"], r["k_sigma"], r["k_G"], r["width_mm"],
                            -(r["sigma_nom_c"] or 0)))
    for i, r in enumerate(out):     # rank within each block
        prev = out[i - 1] if i else None
        same = prev and all(prev[k] == r[k] for k in ("pair", "load", "k_sigma", "k_G", "width_mm", "GIIc_over_GIc"))
        r["rank"] = prev["rank"] + 1 if same else 1
    return out


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--plots", action="store_true")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--materials", default=str(ROOT / "config" / "materials.json"))
    ap.add_argument("--study", default=str(ROOT / "config" / "study.json"))
    ap.add_argument("--out", default=str(ROOT / "results"))
    ap.add_argument("--runs", default=None, help="folder with the setups (default: runs/ of the project)")
    ap.add_argument("--raw", default="raw", help="raw subdirectory name inside each setup")
    ap.add_argument("--benchmark", action="store_true",
                    help="evaluate runs/benchmark_griffith (use --raw raw_synthetic without Abaqus)")
    a = ap.parse_args()

    study = json.loads(pathlib.Path(a.study).read_text(encoding="utf-8"))
    mats = load_materials(a.materials)
    runs = pathlib.Path(a.runs) if a.runs else ROOT / study.get("runs_dir", "runs")
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    plot_dir = out / "curves"
    if a.plots:
        plot_dir.mkdir(exist_ok=True)

    if a.benchmark:
        dirs = [runs / "benchmark_griffith"]
        study = dict(study, post=dict(k_sigma=[1.0], k_G=[1.0], feature_width_mm=[1.0]))
    else:
        dirs = sorted(p.parent for p in runs.glob("*/case.json")
                      if p.parent.name != "benchmark_griffith" and a.only in p.parent.name)
    jobs = [(d, study, mats, a.raw, a.plots, plot_dir) for d in dirs]

    rows, skipped = [], []
    if a.jobs > 1:
        with ProcessPoolExecutor(a.jobs) as ex:
            results = list(ex.map(_safe, jobs))
    else:
        results = [_safe(j) for j in jobs]
    for name, r, msg in results:
        if msg:
            skipped.append((name, msg))
        rows.extend(r)

    if not rows:
        print("no results found (FE not run yet?)")
        for n, m in skipped[:10]:
            print("  skipped", n, "-", m)
        return
    gov = governing_rows(rows)
    write_csv(out / "paths.csv", rows, KEYS_PATH)
    write_csv(out / "governing.csv", gov, KEYS_PATH + ["runner_up", "runner_up_kind", "margin", "complete"])
    rk = ranking_rows(gov)
    write_csv(out / "ranking.csv", rk, list(rk[0].keys()))

    n_set = len({r["setup"] for r in rows})
    flags = [g for g in gov if g["at_lower"] or g["at_upper"]]
    n_lo = sum(1 for g in gov if g["at_lower"]); n_up = sum(1 for g in gov if g["at_upper"])
    n_seg = sum(1 for g in gov if g.get("segment_exhausted"))
    n_cc = sum(r["n_closure_only"] for r in rows if r["k_sigma"] == rows[0]["k_sigma"]
               and r["k_G"] == rows[0]["k_G"] and r["width_mm"] == rows[0]["width_mm"])
    inter = [g for g in gov if g["interpenetration"]]
    badcc = [g for g in gov if np.isfinite(g["cc_ratio_median"]) and not (0.9 < g["cc_ratio_median"] < 1.1)]
    incomplete = sorted({g["setup"] for g in gov if not g["complete"]})
    print(f"{n_set} setups, {len(rows)} path evaluations -> {out}")
    print(f"  governing results with da-bound flag:        {len(flags)}  "
          f"(lower {n_lo}: refine small da / upper {n_up}: longer paths)")
    print(f"  governing: whole segment debonds at onset:   {n_seg}  (physical; next step = sequential analysis)")
    print(f"  crack increments with G from closure only:   {n_cc}  (Delta ALLSE below output precision)")
    eff = [abs(g["lam_c_nopen"] / g["lam_c"] - 1) for g in inter]
    print(f"  governing results with crack-face overlap:   {len(inter)}  (overlap >1% and >1% of G_I energy)"
          + (f"; effect of removing it: median {np.median(eff):.1%}, max {max(eff):.1%}" if eff else ""))
    print(f"  governing results with closure/energy off:   {len(badcc)}  (check mesh / NFORC side)")
    print(f"  setups with missing cracked results:         {len(incomplete)}")
    for n, m in skipped:
        print("  skipped", n, "-", m)


def _safe(job):
    try:
        return process_setup(job)
    except Exception as e:  # keep the sweep going; report at the end
        return job[0].name, [], "ERROR " + repr(e) + "\n" + traceback.format_exc(limit=2)


if __name__ == "__main__":
    main()
