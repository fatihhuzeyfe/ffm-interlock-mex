"""
Verification of the FFM chain against an exact solution.

Griffith crack of half length a in an infinite plate under remote stress sigma,
symmetric growth of both tips by da (b = a + da):
    sigma_yy(x) = sigma * x / sqrt(x^2 - a^2)           (x > a, ahead of the tip)
    G_inc(da)   = pi sigma^2 (a + b) / (2 E')            (per tip)
    [u_y](x)    = 4 sigma / E' * sqrt(b^2 - x^2)         (crack of half length b)
Coupled criterion -> 2 l_ch / (pi (a+b)) = (b^2 - a^2)/b^2, sigma_f^2 = sigma_c^2 (1 - a^2/b^2)
Limits: a -> 0 gives sigma_f -> sigma_c, a -> inf gives Griffith sqrt(E' G_c / (pi a)).

Part 1 feeds the analytical fields to the kernel.
Part 2 writes *synthetic Abaqus raw files* for runs/benchmark_griffith and runs the
full post-processor (traction projection, seam jumps, crack closure, Delta ALLSE,
area normalisation, centred-path folding). If both pass, any discrepancy in the
real FE benchmark is due to the FE model (mesh, partitions, seams), not the code.
"""
import sys, json, math, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffm.kernel import Criterion, ffm_path

E_p, sigma_c, G_c = 3000.0, 50.0, 3.0            # any consistent values
l_ch = E_p * G_c / sigma_c ** 2


def exact(a):
    f = lambda b: 2 * l_ch / (math.pi * (a + b)) - (b * b - a * a) / (b * b)
    b = brentq(f, a * (1 + 1e-12) + 1e-15, a + 1e3 * l_ch + 10 * a)
    return sigma_c * math.sqrt(1 - a * a / (b * b)), b - a


def kernel_solution(a, n_da=14, n_r=4000):
    crit = Criterion(sigma_c=sigma_c, G_Ic=G_c)
    da_hi = 20 * l_ch
    r = np.concatenate([[0.0], np.geomspace(1e-6 * max(a, l_ch), da_hi * 1.01, n_r)])
    x = a + r
    with np.errstate(divide="ignore"):
        syy = np.where(r > 0, x / np.sqrt(np.maximum(x * x - a * a, 1e-300)), 1e30)
    da = np.geomspace(1e-4 * l_ch, da_hi, n_da)
    GI = math.pi * (2 * a + da) / (2 * E_p)
    res = ffm_path([(r, syy, 0 * r)], da, GI, 0 * GI, crit)
    return res.lam_c, res.da_c, res


def part1(plot=True):
    a_list = np.geomspace(1e-2, 1e2, 25) * l_ch
    rows = []
    for a in a_list:
        sf, dac = exact(a)
        lam, da_k, _ = kernel_solution(a)
        rows.append((a / l_ch, sf / sigma_c, lam / sigma_c, dac / l_ch, da_k / l_ch))
    rows = np.array(rows)
    err_s = np.abs(rows[:, 2] / rows[:, 1] - 1)
    err_a = np.abs(rows[:, 4] / rows[:, 3] - 1)
    print(f"[part 1] kernel vs exact: max |error| strength {err_s.max():.2e}, onset length {err_a.max():.2e}")
    if plot:
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
        aa = np.geomspace(1e-2, 1e2, 200)
        ax[0].loglog(aa, [exact(x * l_ch)[0] / sigma_c for x in aa], "k-", label="exact FFM")
        ax[0].loglog(rows[:, 0], rows[:, 2], "o", mfc="none", label="kernel")
        ax[0].loglog(aa, np.sqrt(1 / (math.pi * aa)), "--", color="gray", label="Griffith (LEFM)")
        ax[0].axhline(1, ls=":", color="gray"); ax[0].set_ylim(0.03, 1.5)
        ax[0].set_xlabel("a / l_ch"); ax[0].set_ylabel("sigma_f / sigma_c"); ax[0].legend(fontsize=8)
        _, _, r = kernel_solution(0.5 * l_ch)
        r2 = ffm_path(*_curves_args(0.5 * l_ch), return_curves=True)
        c = r2.curves
        ax[1].semilogx(c["da"] / l_ch, c["lam_s"] / sigma_c, label="stress condition")
        ax[1].semilogx(c["da"] / l_ch, c["lam_G"] / sigma_c, label="energy condition")
        ax[1].plot(r2.da_c / l_ch, r2.lam_c / sigma_c, "ko", label="onset (min of max)")
        ax[1].set_ylim(0, 1.2); ax[1].set_xlabel("da / l_ch"); ax[1].set_ylabel("load / sigma_c")
        ax[1].set_title("a = 0.5 l_ch", fontsize=9); ax[1].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(ROOT / "docs/benchmark_griffith_kernel.png", dpi=140)
    return err_s.max(), err_a.max()


def _curves_args(a):
    crit = Criterion(sigma_c=sigma_c, G_Ic=G_c)
    r = np.concatenate([[0.0], np.geomspace(1e-6 * l_ch, 20 * l_ch, 4000)])
    x = a + r
    with np.errstate(divide="ignore"):
        syy = np.where(r > 0, x / np.sqrt(np.maximum(x * x - a * a, 1e-300)), 1e30)
    da = np.geomspace(1e-4 * l_ch, 20 * l_ch, 14)
    GI = math.pi * (2 * a + da) / (2 * E_p)
    return [(r, syy, 0 * r)], da, GI, 0 * GI, crit


# ---------------------------------------------------------------------------
def part2():
    """Synthetic raw files -> full post-processor."""
    from ffm.postprocess import Setup, evaluate
    d = ROOT / "runs" / "benchmark_griffith"
    case = json.loads((d / "case.json").read_text())
    a = case["benchmark"]["a"]
    W = case["W"]; xc = W / 2
    sig_ref = 1.0
    raw = d / "raw_synthetic"
    raw.mkdir(exist_ok=True)
    # "mesh" nodes on the interface: geometric grading around the tips, like the FE seeds
    x_nodes = np.unique(np.concatenate([
        xc + a + np.concatenate([[0], np.geomspace(1e-5, W / 2 - a - 1e-3, 1500)]),
        xc - a - np.concatenate([[0], np.geomspace(1e-5, W / 2 - a - 1e-3, 1500)]),
        np.linspace(xc - a, xc + a, 200)]))
    # lumped nodal forces across the line: F_i = int sigma_yy over the tributary length
    #   int sigma r / sqrt(r^2 - a^2) dr = sigma sqrt(r^2 - a^2)  (r > a), 0 on the pre-crack
    prim = lambda r: sig_ref * math.sqrt(max(r * r - a * a, 0.0)) if r > a else 0.0
    mids = np.concatenate([[x_nodes[0]], 0.5 * (x_nodes[1:] + x_nodes[:-1]), [x_nodes[-1]]])
    rows = []
    for i, x in enumerate(x_nodes):
        r = abs(x - xc)
        syy = sig_ref * r / math.sqrt(r * r - a * a) if r > a + 1e-12 else (1e6 if abs(r - a) < 1e-12 else 0.0)
        r_lo, r_hi = sorted([abs(mids[i] - xc), abs(mids[i + 1] - xc)])
        if (mids[i] - xc) * (mids[i + 1] - xc) < 0:      # tributary interval straddles the centre
            Fy = 0.0
        else:
            Fy = prim(r_hi) - prim(r_lo)
        s = [0.0, syy, 0.0]
        rows.append([x, 0, x, 0.0, s, s, [0.0, Fy]])
    U0 = sig_ref ** 2 / (2 * E_p) * W * 2 * case["H"]      # stored energy of the uniformly loaded plate
    (raw / "uncracked.json").write_text(json.dumps(dict(ALLSE=U0, ALLWK=U0, F=sig_ref * W, interface=rows, bulk={})))
    path = case["paths"][0]
    for k, cr in enumerate(path["cracks"]):
        b = a + cr["da"]
        xs = x_nodes[(np.abs(x_nodes - xc) < b - 1e-12)]
        jumps = [[x, 0, 0.0, 4 * sig_ref / E_p * math.sqrt(b * b - (x - xc) ** 2)] for x in xs]
        dU = math.pi * sig_ref ** 2 * (b * b - a * a) / E_p          # both tips, per unit thickness
        (raw / f"{path['pid']}__{k:02d}.json").write_text(json.dumps(dict(ALLSE=U0 - dU, ALLWK=0, F=0, jumps=jumps)))
    # run the post-processor on the synthetic data
    st = Setup(d, raw_subdir="raw_synthetic")
    pd = st.path_data(path)
    cc = pd.cc_ratio
    print(f"[part 2] crack-closure / energy ratio: min {cc.min():.4f}, max {cc.max():.4f}")
    crit = Criterion(sigma_c=sigma_c, G_Ic=G_c)
    worst = 0.0
    # physical size scaling a_phys = kappa * a. The da grid of the case (0.01a..2a) only
    # brackets the onset length for a/l_ch >~ 0.2; below that the optimum sits on the upper
    # bound and must be flagged (this is the check that protects the parametric study).
    for kappa in (0.1, 1.0, 3.0, 10.0, 30.0):
        r0 = evaluate(st, crit, {}, kappa=kappa, paths=[pd])[0]
        sf_num = r0["lam_c"] * sig_ref
        sf_ex, dac = exact(a * kappa)
        flag = "UPPER-BOUND FLAG" if r0["at_upper"] else ("LOWER-BOUND FLAG" if r0["at_lower"] else "")
        if not (r0["at_upper"] or r0["at_lower"]):
            worst = max(worst, abs(sf_num / sf_ex - 1))
        print(f"   a/l_ch = {a * kappa / l_ch:6.3f}: FFM chain {sf_num:7.3f} MPa, exact {sf_ex:7.3f} MPa, "
              f"da_c {r0['da_c']:.4f} vs {dac:.4f} mm  {flag}")
    print(f"[part 2] max |error| {worst:.2e}")
    return cc, worst


if __name__ == "__main__":
    e1 = part1()
    cc, e2 = part2()
    ok = e1[0] < 5e-3 and e2 < 1e-2 and abs(np.median(cc) - 1) < 1e-2
    print("BENCHMARK", "PASSED" if ok else "FAILED")
