"""Preview interface geometries and candidate crack paths (sanity check before FE)."""
import sys, pathlib, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ffm.geometry import (GeomSpec, interface_polyline, candidate_paths, mirror_prune,
                          default_da, arc_length, point_at, singular_points, sine_centres)


import json
_study = json.loads((pathlib.Path(__file__).resolve().parents[1] / "config" / "study.json").read_text())
_da = _study["da"]
EXT = dict(n=_da["n"], min=_da["min"], max=_da["max_ext"]) if _da.get("extend") else None


def draw(ax, g, prune=True):
    P = interface_polyline(g); S = arc_length(P)
    paths = candidate_paths(g, P, default_da(), centred_S=sine_centres(g, P) if g.shape == "sine" else None,
                            bulk_offsets=_study["bulk_offsets_deg"], extend=EXT,
                            periodic=_study.get("lateral") == "periodic")
    if prune:
        paths = mirror_prune(paths, P, g.W)
    ax.fill(list(P[:, 0]) + [g.W, 0], list(P[:, 1]) + [-1.2, -1.2], color="#dbe7f3", zorder=0)
    ax.fill(list(P[:, 0]) + [g.W, 0], list(P[:, 1]) + [3, 3], color="#f6e3cf", zorder=0)
    ax.plot(P[:, 0], P[:, 1], "k-", lw=1.2)
    for p in paths:
        if p.kind == "interface":
            a, b = p.crack_interval(p.L)
            xy = np.array([point_at(P, S, s) for s in np.linspace(a, b, 40)])
            ax.plot(xy[:, 0], xy[:, 1], "-", color="#c0392b", lw=3, alpha=.55)
            if p.mode == "centred":
                ax.plot(*point_at(P, S, p.S0), "s", ms=4, color="#c0392b")
        else:
            o = p.origin; e = (math.cos(math.radians(p.angle_deg)), math.sin(math.radians(p.angle_deg)))
            ax.plot([o[0], o[0] + p.L * e[0]], [o[1], o[1] + p.L * e[1]], "-", color="#1f6f8b", lw=1.5)
    for sp in singular_points(P):
        ax.plot(*sp.xy, "o", ms=4, color="k")
    ax.set_title(f"{g.tag()}  ({len(paths)} paths)", fontsize=9)
    ax.set_aspect("equal"); ax.set_xlim(-0.2, g.W + 0.2); ax.set_ylim(-1.2, 3.0); ax.axis("off")


if __name__ == "__main__":
    specs = [GeomSpec("flat"), GeomSpec("U", ar=1.0), GeomSpec("T", ar=1.0),
             GeomSpec("dovetail", ar=1.0, alpha_deg=20), GeomSpec("dovetail", ar=2.0, alpha_deg=30),
             GeomSpec("sine", ar=1.0)]
    fig, axs = plt.subplots(2, 3, figsize=(12, 6))
    for ax, g in zip(axs.flat, specs):
        draw(ax, g)
    fig.suptitle("A (blue, below) / B (orange, above); red = interface paths, teal = bulk paths "
                 "(mirror-pruned, one symmetric half shown)", fontsize=10)
    out = sys.argv[1] if len(sys.argv) > 1 else "docs/geometries.png"
    fig.tight_layout(); fig.savefig(out, dpi=140); print("saved", out)
