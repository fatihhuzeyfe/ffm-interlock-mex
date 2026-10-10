"""Mesh convergence study for the paper (verification D3).

Three systematically refined meshes with constant refinement ratio r = 2 near the crack lines
(h_line and the number of elements per crack increment are refined together; the far-field
size h_global is kept). The quantities of interest are the FFM onset stress of the governing
mechanism (and, optionally, of every candidate path) for all k_sigma x k_G x w x G_IIc/G_Ic
combinations. Discretisation uncertainty is reported with the Grid Convergence Index
(Roache; procedure of Celik et al., J. Fluids Eng. 130 (2008) 078001), plus two-level
changes and mechanism agreement.

1) prepare the meshes (writes config/study_mesh_<label>.json and prints the commands)

    python scripts/mesh_convergence.py prepare            # representative set (see REPRESENTATIVE)

   Levels (default): coarse  h_line 0.10, n_per_edge 2  -> runs_mesh_coarse / results_mesh_coarse
                     medium  h_line 0.05, n_per_edge 4  -> runs / results            (production mesh)
                     fine    h_line 0.025, n_per_edge 8 -> runs_fine / results_fine
   Existing medium/fine results are reused; only the missing level is computed.

2) analyse (after all levels are post-processed)

    python scripts/mesh_convergence.py analyze

   Writes docs/mesh_convergence/: convergence_table.csv, summary.md (paper-ready numbers and a
   LaTeX table), fig_mesh_convergence.png/.pdf.
"""
import sys, json, csv, math, argparse, pathlib, statistics as st, re

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

LEVELS = [  # label, h_line, n_per_edge, runs dir, results dir   (coarse -> fine)
    ("coarse", 0.10, 2, "runs_mesh_coarse", "results_mesh_coarse"),
    ("medium", 0.05, 4, "runs", "results"),
    ("fine", 0.025, 8, "runs_fine", "results_fine"),
]
KEY = ("setup", "k_sigma", "k_G", "width_mm", "GIIc_over_GIc")
FS = 1.25            # safety factor for three-grid studies (Roache / Celik et al.)


# ----------------------------------------------------------------------------- prepare
REPRESENTATIVE = [                 # one geometry per family, all three pairs, both loads
    "dovetail_ar1_a20__PLA-PETG",       # tension + shear, r0 + r90 (fine level already computed)
    "T_ar1__PLA-ABS__tension__r0",      # several re-entrant corners, crossing bulk lines
    "dovetail_ar2_a10__PLA-ABS__tension__r0",
    "U_ar1__PETG-ABS__tension__r0",     # smooth arc + corners
    "sine_ar1__PETG-ABS__shear__r0",    # smooth interface, centred paths, mode II
    "flat__PLA-ABS__tension__r0",       # nucleation in a uniform field
]


def prepare(a):
    base = json.loads((ROOT / "config" / "study.json").read_text(encoding="utf-8"))
    sel = [x for x in a.setups.split(",") if x] if a.setups else REPRESENTATIVE
    print("Medium level = production runs/results (from the main sweep).")
    print("Run these commands for the levels that are not complete yet:\n")
    for label, hl, n, runs, res in LEVELS:
        s = json.loads(json.dumps(base))
        s["runs_dir"] = runs
        s["mesh"] = dict(s["mesh"], h_line=hl, n_per_edge=n)
        s["_notes"] = f"mesh convergence level '{label}' (only the mesh differs from study.json)"
        f = ROOT / "config" / f"study_mesh_{label}.json"
        f.write_text(json.dumps(s, indent=1), encoding="utf-8")
        if runs == "runs":
            continue                                  # production level: comes from the sweep
        print(f"# ---- {label}: h_line = {hl}, n_per_edge = {n}")
        for x in sel:
            print(f"python scripts/make_cases.py --study config/study_mesh_{label}.json --only {x}")
        for x in sel:
            print(f"python scripts/launch_abaqus.py --runs {runs} --only {x} --parallel {a.parallel}")
        print(f"python scripts/run_postprocess.py --study config/study_mesh_{label}.json --out {res}\n")
    print("# ---- after the production post-processing (results/governing.csv) exists:")
    print("python scripts/mesh_convergence.py analyze")


# ----------------------------------------------------------------------------- analysis helpers
def load_gov(res):
    f = ROOT / res / "governing.csv"
    if not f.exists():
        return None
    return {tuple(r[k] for k in KEY): r for r in csv.DictReader(open(f, encoding="utf-8"))}


_case_cache = {}


def mechanism(runs, setup, pid):
    """Mechanism key independent of the path label: interface edge paths are identified by the
    arc-length interval they cover (I2f and I3b on the same head-top segment are the same
    mechanism), bulk paths by origin and angle."""
    k = (runs, setup)
    if k not in _case_cache:
        f = ROOT / runs / setup / "case.json"
        _case_cache[k] = {p["pid"]: p for p in json.loads(f.read_text())["paths"]} if f.exists() else {}
    p = _case_cache[k].get(pid)
    if p is None:
        return pid
    if p["kind"] == "interface" and p["mode"] == "edge":
        s0, s1 = sorted((p["S0"], p["S0"] + p["direction"] * p["L"]))
        return f"interface[{s0:.3f},{s1:.3f}]"
    if p["kind"] == "interface":
        return f"interface@{p['S0']:.3f}"
    return f"bulk {p['material']} {p['origin'][0]:.3f},{p['origin'][1]:.3f} {p['angle_deg']:.0f}"


def element_counts(runs, setups):
    n = []
    for s in setups:
        f = ROOT / runs / s / "abaqus_log.txt"
        if f.exists():
            m = re.findall(r"uncracked: (\d+) elements", f.read_text(errors="ignore"))
            if m:
                n.append(int(m[-1]))
    return st.mean(n) if n else float("nan")


def gci3(f1, f2, f3, r=2.0):
    """f1 fine, f2 medium, f3 coarse. Returns dict with observed order p, extrapolated value,
    GCI_fine (relative), convergence type. Constant refinement ratio r (Celik et al. 2008)."""
    e21, e32 = f2 - f1, f3 - f2
    scale = max(abs(f1), 1e-300)
    if abs(e21) < 1e-6 * scale and abs(e32) < 1e-6 * scale:
        return dict(p=float("nan"), f_ext=f1, gci_fine=0.0, kind="converged (no change)")
    if abs(e21) < 1e-6 * scale:
        return dict(p=float("nan"), f_ext=f1, gci_fine=0.0, kind="fine = medium")
    ratio = e32 / e21
    kind = "monotonic" if ratio > 0 else "oscillatory"
    p = abs(math.log(abs(ratio))) / math.log(r) if abs(ratio) != 1 else float("nan")
    if not math.isfinite(p) or p <= 0.05:
        # no asymptotic behaviour: conservative estimate with p = 1
        p_use, kind = 1.0, kind + " (p<=0: p=1 used)"
    else:
        p_use = min(p, 4.0)                      # cap unrealistic orders
    f_ext = (r ** p_use * f1 - f2) / (r ** p_use - 1)
    gci = FS * abs(e21 / f1) / (r ** p_use - 1)
    return dict(p=p, f_ext=f_ext, gci_fine=gci, kind=kind)


def pct(v):
    return f"{100 * v:.2f} %"


# ----------------------------------------------------------------------------- analyze
def analyze(a):
    levels = [(lab, hl, n, runs, res, load_gov(res)) for lab, hl, n, runs, res in LEVELS]
    have = [L for L in levels if L[5] is not None]
    if len(have) < 2:
        sys.exit("need at least two post-processed levels (see 'prepare')")
    out = ROOT / "docs" / "mesh_convergence"
    out.mkdir(parents=True, exist_ok=True)
    common = set.intersection(*[set(L[5]) for L in have])
    if a.only:
        flt = a.only.split(",")
        common = {k for k in common if any(x in k[0] for x in flt)}
    common = sorted(common)
    setups = sorted({k[0] for k in common})
    print(f"levels: {[L[0] for L in have]}, {len(common)} common results, {len(setups)} setups")

    rows = []
    fine, med = have[-1], have[-2]
    coarse = have[-3] if len(have) >= 3 else None
    for k in common:
        f1 = float(fine[5][k]["sigma_nom_c"]); f2 = float(med[5][k]["sigma_nom_c"])
        mech = [mechanism(L[3], k[0], L[5][k]["pid"]) for L in have]
        r = dict(zip(KEY, k), sigma_fine=f1, sigma_medium=f2,
                 change_medium_fine=f2 / f1 - 1, mechanism_fine=mech[-1],
                 same_mechanism=all(m == mech[-1] for m in mech),
                 same_kind=len({L[5][k]["kind"] for L in have}) == 1)
        if coarse:
            f3 = float(coarse[5][k]["sigma_nom_c"])
            g = gci3(f1, f2, f3)
            r.update(sigma_coarse=f3, change_coarse_fine=f3 / f1 - 1, p_observed=g["p"],
                     sigma_extrapolated=g["f_ext"], gci_fine=g["gci_fine"],
                     gci_medium=g["gci_fine"] * 2 ** (g["p"] if math.isfinite(g["p"]) else 1.0),
                     convergence=g["kind"])
        rows.append(r)

    keys = list(rows[0].keys())
    with open(out / "convergence_table.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)

    # ---------------- summary numbers
    ch = [abs(r["change_medium_fine"]) for r in rows]
    n_el = {L[0]: element_counts(L[3], setups) for L in have}
    lines = ["# Mesh convergence (verification D3)\n",
             f"Setups: {', '.join(setups)}  \nResults compared per level: {len(rows)} "
             f"(k_sigma x k_G x w x G_IIc/G_Ic combinations)\n",
             "| Level | h_line [mm, w = 1] | elements per crack increment | mean elements (uncracked) |",
             "|---|---|---|---|"]
    for L in have:
        lines.append(f"| {L[0]} | {L[1]} | {L[2]} | {n_el[L[0]]:.0f} |")
    lines += ["", "## Production (medium) mesh vs fine mesh\n",
              f"- onset stress change: median {pct(st.median([r['change_medium_fine'] for r in rows]))}, "
              f"mean |.| {pct(st.mean(ch))}, 95th percentile {pct(sorted(ch)[int(0.95 * (len(ch) - 1))])}, "
              f"max {pct(max(ch))}",
              f"- same failure mechanism (all levels): {sum(r['same_mechanism'] for r in rows)}/{len(rows)}; "
              f"same kind (interface/bulk): {sum(r['same_kind'] for r in rows)}/{len(rows)}"]
    if coarse:
        g = [r["gci_fine"] for r in rows]
        gm = [r["gci_medium"] for r in rows]
        ps = [r["p_observed"] for r in rows if math.isfinite(r["p_observed"])]
        mono = sum(r["convergence"].startswith("monotonic") for r in rows)
        lines += ["", "## Three-grid analysis (GCI, Fs = 1.25, r = 2)\n",
                  f"- GCI of the production mesh: median {pct(st.median(gm))}, 95th percentile "
                  f"{pct(sorted(gm)[int(0.95 * (len(gm) - 1))])}, max {pct(max(gm))}",
                  f"- GCI of the fine mesh: median {pct(st.median(g))}, max {pct(max(g))}",
                  f"- observed order of convergence p: median {st.median(ps):.2f} (n = {len(ps)})" if ps else
                  "- observed order: not defined (changes below 1e-6)",
                  f"- monotonic convergence: {mono}/{len(rows)}; others oscillatory or already converged"]
    # per setup
    lines += ["", "## Per setup\n", "| Setup | mean abs change medium-fine | max | same mechanism |",
              "|---|---|---|---|"]
    for s in setups:
        rs = [r for r in rows if r["setup"] == s]
        c = [abs(r["change_medium_fine"]) for r in rs]
        lines.append(f"| {s} | {pct(st.mean(c))} | {pct(max(c))} | "
                     f"{sum(r['same_mechanism'] for r in rs)}/{len(rs)} |")
    # LaTeX
    lines += ["", "## LaTeX table\n", "```latex",
              r"\begin{tabular}{lccc}", r"\hline",
              r"Mesh & $h_\mathrm{line}/w$ & elements per $\Delta a$ & elements \\ \hline"]
    for L in have:
        lines.append(f"{L[0]} & {L[1]} & {L[2]} & {n_el[L[0]]:.0f} \\\\")
    lines += [r"\hline", r"\end{tabular}", "```", "",
              "Suggested text: *Refining the mesh along all candidate crack paths by a factor of two "
              f"changed the predicted onset stress by {pct(st.mean(ch))} on average (maximum {pct(max(ch))}); "
              f"the predicted onset mechanism was unchanged in {sum(r['same_mechanism'] for r in rows)} of "
              f"{len(rows)} cases.*" + (f" *The grid convergence index of the production mesh is "
                                        f"{pct(st.median([r['gci_medium'] for r in rows]))} (median).*" if coarse else "")]
    (out / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[:20]))

    # ---------------- figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    hs = [L[1] for L in have]
    for s in setups:
        rs = [r for r in rows if r["setup"] == s and float(r["k_sigma"]) == 1.0 and float(r["k_G"]) == 1.0
              and float(r["GIIc_over_GIc"]) == 2.0]
        for r in rs:
            vals = ([r["sigma_coarse"]] if coarse else []) + [r["sigma_medium"], r["sigma_fine"]]
            ref = r.get("sigma_extrapolated", r["sigma_fine"])
            lbl = s.replace("dovetail_", "").replace("__", " ") + f", w={float(r['width_mm']):g}"
            ax[0].plot(hs, [v / ref for v in vals], "o-", ms=4, lw=1, label=lbl)
    ax[0].set_xscale("log"); ax[0].invert_xaxis()
    ax[0].set_xlabel("element size along crack lines $h_{line}/w$")
    ax[0].set_ylabel(r"$\sigma_c / \sigma_{c,ext}$" if coarse else r"$\sigma_c / \sigma_{c,fine}$")
    ax[0].axhline(1, color="k", lw=0.5)
    ax[0].legend(fontsize=6, ncol=2)
    ax[0].set_title("(a) onset stress vs mesh (k = 1)", fontsize=9)
    data = [[100 * abs(r["change_medium_fine"]) for r in rows if r["setup"] == s] for s in setups]
    ax[1].boxplot(data, vert=True)
    ax[1].set_xticks(range(1, len(setups) + 1))
    ax[1].set_xticklabels([s.split("__", 1)[1].replace("__", "\n") for s in setups], fontsize=7)
    ax[1].set_ylabel("|change| production vs fine mesh [%]")
    ax[1].set_title(f"(b) all {len(rows)} parameter combinations", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fig_mesh_convergence.{ext}", dpi=200)
    print(f"\nwritten: {out}")


def make(a):
    """Generate the coarse and fine case files of the representative set (no Abaqus)."""
    import subprocess
    base = json.loads((ROOT / "config" / "study.json").read_text(encoding="utf-8"))
    sel = [x for x in a.setups.split(",") if x] if a.setups else REPRESENTATIVE
    for label, hl, n, runs, res in LEVELS:
        s = json.loads(json.dumps(base))
        s["runs_dir"] = runs
        s["mesh"] = dict(s["mesh"], h_line=hl, n_per_edge=n)
        s["_notes"] = f"mesh convergence level '{label}' (only the mesh differs from study.json)"
        f = ROOT / "config" / f"study_mesh_{label}.json"
        f.write_text(json.dumps(s, indent=1), encoding="utf-8")
        if runs == "runs":
            continue
        for x in sel:
            subprocess.run([sys.executable, str(ROOT / "scripts" / "make_cases.py"), "--study", str(f),
                            "--only", x], check=True, stdout=subprocess.DEVNULL)
        print(f"{label}: cases written to {runs}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("prepare"); p1.add_argument("--setups", default="",
                                                    help="comma-separated setup filters (default: representative set)")
    p1.add_argument("--parallel", type=int, default=4)
    p2 = sub.add_parser("analyze"); p2.add_argument("--only", default="")
    p3 = sub.add_parser("make"); p3.add_argument("--setups", default="")
    a = ap.parse_args()
    {"prepare": prepare, "analyze": analyze, "make": make}[a.cmd](a)


if __name__ == "__main__":
    main()
