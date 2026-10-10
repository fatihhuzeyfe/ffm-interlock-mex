"""Compare two result folders (e.g. mesh convergence: results vs results_fine).

    python scripts/compare_results.py results results_fine

Matches governing rows on setup x k_sigma x k_G x width x G_IIc/G_Ic and reports the relative
change of the onset stress and whether the governing path / mechanism stays the same.
Mesh is converged when |change| < 2 % (median) and the governing mechanism agrees.
"""
import sys, csv, statistics as st

def load(d):
    rows = list(csv.DictReader(open(f"{d}/governing.csv", encoding="utf-8")))
    key = lambda r: (r["setup"], r["k_sigma"], r["k_G"], r["width_mm"], r["GIIc_over_GIc"])
    return {key(r): r for r in rows}

a, b = load(sys.argv[1]), load(sys.argv[2])
common = sorted(set(a) & set(b))
if not common:
    sys.exit("no common rows")
d = [float(b[k]["sigma_nom_c"]) / float(a[k]["sigma_nom_c"]) - 1 for k in common]
same_pid = sum(a[k]["pid"] == b[k]["pid"] for k in common)
same_kind = sum(a[k]["kind"] == b[k]["kind"] for k in common)
print(f"{len(common)} common results")
print(f"onset stress change: median {st.median(d):+.2%}, mean |.| {st.mean(map(abs, d)):.2%}, "
      f"max |.| {max(map(abs, d)):.2%}")
print(f"same governing path: {same_pid}/{len(common)}, same kind (interface/bulk): {same_kind}/{len(common)}")
by = {}
for k, x in zip(common, d):
    by.setdefault(k[0], []).append(abs(x))
for s, v in sorted(by.items()):
    print(f"  {s:50s} mean |change| {st.mean(v):.2%}")
