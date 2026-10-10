"""Publication contour plots of FE fields (stress, opening) without Abaqus/Viewer.

    python scripts/fe_figure.py <setup>                         # uncracked model
    python scripts/fe_figure.py <setup> --crack I2f__11         # a cracked configuration
           [--field S22|S11|S12|s1|smax_ratio] [--window x0 x1 y0 y1] [--scale 50]
           [--abaqus abaqus] [--out docs/figures]

Steps (each skipped if already done):
  1. the ODB: uncracked -> runs/<setup>/work_s0/uncr.odb (kept by every sweep run);
     cracked  -> re-run of that single analysis with its ODB kept (driver option keep=...)
  2. export: abaqus python abaqus/export_field.py -> runs/<setup>/fields/<name>.json
  3. plot  : matplotlib, both materials with their own nodal averages (stress jump at the
     interface preserved), interface and crack faces drawn, optional deformed shape.
Fields: S11, S22, S12 (global axes, MPa per reference load), s1 = max principal stress.
Values are for the reference displacement of the FE model; multiply by lambda_c (paths.csv)
to show the field at onset.
"""
import sys, json, argparse, pathlib, subprocess, os
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]


def sh(cmd):
    print("$", cmd)
    return subprocess.call(cmd, shell=True, cwd=str(ROOT))


def ensure_field(setup, crack, abq):
    d = ROOT / "runs" / setup
    tag = crack or "uncracked"
    fj = d / "fields" / f"{tag}.json"
    if fj.exists():
        return fj
    work = d / "work_s0"
    odb = work / ("uncr.odb" if not crack else f"c_{crack.replace('.', '_')}.odb")
    if not odb.exists():
        drv = ROOT / "abaqus" / "build_and_run.py"
        extra = f"keep={crack}" if crack else "keep=uncracked"
        sh(f'{abq} cae noGUI="{drv}" -- "{d / "case.json"}" 0 1 {extra}')
    if not odb.exists():
        sys.exit(f"ODB not found: {odb}")
    fj.parent.mkdir(exist_ok=True)
    sh(f'{abq} python "{ROOT / "abaqus" / "export_field.py"}" "{odb}" "{d / "case.json"}" "{fj}"')
    if not fj.exists():
        sys.exit("export failed")
    return fj


def split_tri6(t):
    """Abaqus 6-node triangle (1,2,3 corners; 4,5,6 mid-sides 12,23,31) -> 4 linear triangles."""
    a, b, c, ab, bc, ca = t
    return [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]


def field_values(st, name):
    s11, s22, s12 = st[:, 0], st[:, 1], st[:, 2]
    if name == "S11":
        return s11
    if name == "S22":
        return s22
    if name == "S12":
        return s12
    if name == "s1":
        return 0.5 * (s11 + s22) + np.sqrt((0.5 * (s11 - s22)) ** 2 + s12 ** 2)
    raise SystemExit(f"unknown field {name}")


def plot(fj, case, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    data = json.loads(fj.read_text())
    fig, ax = plt.subplots(figsize=(6.5, 6.5 * (a.window[3] - a.window[2]) / (a.window[1] - a.window[0]))
                           if a.window else (5, 8))
    vals_all = []
    tris = []
    for key, m in data["materials"].items():
        xy = np.array(m["xy"]); u = np.array(m["u"]); st = np.array(m["stress"])
        xy = xy + a.scale * u
        tri = np.array([s for t in m["tri6"] for s in split_tri6(t)])
        v = field_values(st, a.field)
        tris.append((mtri.Triangulation(xy[:, 0], xy[:, 1], tri), v))
        vals_all.append(v if not a.window else v[(xy[:, 0] >= a.window[0]) & (xy[:, 0] <= a.window[1]) &
                                                (xy[:, 1] >= a.window[2]) & (xy[:, 1] <= a.window[3])])
    allv = np.concatenate([x for x in vals_all if len(x)])
    lo, hi = np.percentile(allv, [1, 99])
    if a.field in ("S12",) or lo < 0 < hi:
        m_ = max(abs(lo), abs(hi)); lo, hi = -m_, m_
        cmap = "RdBu_r"
    else:
        cmap = "viridis"
    levels = np.linspace(lo, hi, 41)
    for T, v in tris:
        cs = ax.tricontourf(T, np.clip(v, lo, hi), levels=levels, cmap=cmap)
    P = np.array(case["interface_points"])
    ax.plot(P[:, 0], P[:, 1], "k-", lw=0.6)
    ax.set_aspect("equal")
    if a.window:
        ax.set_xlim(a.window[0], a.window[1]); ax.set_ylim(a.window[2], a.window[3])
    cb = fig.colorbar(cs, ax=ax, shrink=0.8)
    cb.set_label(f"{a.field} [MPa] (reference load)")
    ttl = case["name"] + ("" if not a.crack else f"  crack {a.crack}")
    ax.set_title(ttl, fontsize=8)
    ax.set_xlabel("x / w"); ax.set_ylabel("y / w")
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    stem = f"{case['name']}__{a.crack or 'uncracked'}__{a.field}"
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    print("figure:", out / f"{stem}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("setup")
    ap.add_argument("--crack", default="")
    ap.add_argument("--field", default="S22")
    ap.add_argument("--window", type=float, nargs=4, default=None)
    ap.add_argument("--scale", type=float, default=0.0, help="deformation scale factor")
    ap.add_argument("--abaqus", default="abaqus")
    ap.add_argument("--out", default=str(ROOT / "docs" / "figures"))
    a = ap.parse_args()
    case = json.loads((ROOT / "runs" / a.setup / "case.json").read_text())
    fj = ensure_field(a.setup, a.crack, a.abaqus)
    plot(fj, case, a)


if __name__ == "__main__":
    main()
