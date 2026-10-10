"""Evaluate the cohesive-zone runs (runs_czm) and compare them with the FFM predictions.

    python scripts/czm_postprocess.py [--runs runs_czm] [--out results_czm] [--fig figures/data]

Per run (nominal stress sigma = F / W of the unit cell):
  sigma_onset          CZM onset = first point where the response deviates by 5 % from the
                       initial linear response (secant stiffness 5 % below the initial one;
                       the non-linearity criterion of fracture testing, cf. ASTM D5528). It
                       captures crack formation both when the load drops (unstable debonding)
                       and when it does not (load transferred to mechanical locking)
  sigma_first_failure  stress when the first interface point is fully damaged
  sigma_instability    first local load maximum followed by a drop of more than DROP (unstable
                       crack growth) = CZM counterpart of the FFM onset
  sigma_separation     stress when the whole interface has failed (>= 98 % of the cohesive nodes
                       fully damaged); beyond it only mechanical locking of the elastic bulk
                       carries load, whose limit (bulk fracture) is not modelled
  reserve              sigma_separation / sigma_instability = load carried by the partially
                       debonded interface after the first instability
  sigma_max            maximum load reached in the run (for information)
Outputs: results_czm/czm_summary.csv; figure data fig5a_czm_curves.csv, fig5b_czm_summary.csv.
"""
import sys, json, argparse, pathlib
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
DROP = 0.02


def analyse(d):
    r = json.loads((d / "raw" / "czm.json").read_text())
    case = json.loads((d / "case.json").read_text())
    cz = case["czm"]
    W = r["W"]
    fr = pd.DataFrame(r["frames"])
    if fr.empty:
        return None, None
    fr["sigma"] = np.abs(fr.F) / W
    fr["u"] = np.abs(fr.u)
    s = fr.sigma.values
    # first instability: first local maximum followed by a drop > DROP before the load recovers
    inst = None
    for i in range(1, len(s) - 1):
        if s[i] >= s[i - 1] and s[i] > s[i + 1]:
            j = i + 1
            while j < len(s) and s[j] <= s[i]:
                if s[j] < (1 - DROP) * s[i]:
                    inst = i
                    break
                j += 1
            if inst is not None:
                break
    i_max = int(np.argmax(s))
    # 5 % deviation from the initial linear response
    u = fr.u.values
    early = (u > 0) & (u <= 0.3 * cz["ffm"]["delta_onset"])
    k0 = float(np.median(s[early] / u[early])) if early.any() else s[1] / max(u[1], 1e-30)
    dev = np.where((u > 0) & (s < 0.95 * k0 * u))[0]
    i5 = int(dev[0]) if len(dev) else i_max
    if i5 > 0:                                   # interpolate the crossing of the 5 % line
        f0 = s[i5 - 1] - 0.95 * k0 * u[i5 - 1]; f1 = s[i5] - 0.95 * k0 * u[i5]
        t = f0 / (f0 - f1) if f0 != f1 else 1.0
        sig5 = float(s[i5 - 1] + t * (s[i5] - s[i5 - 1]))
    else:
        sig5 = float(s[i5])
    if inst is None:
        inst = i_max
    fcol = fr["frac_failed"] if "frac_failed" in fr else pd.Series(np.nan, index=fr.index)
    ff = fr.index[fcol > 0]
    first_fail = float(fr.sigma[ff[0]]) if len(ff) else np.nan
    sep = fr.index[fcol >= 0.98]
    # the last load maximum before (or at) complete separation
    sig_sep = float(fr.sigma[:sep[0] + 1].max()) if len(sep) else np.nan
    ffm = cz["ffm"]
    row = dict(name=case["name"], setup=cz["source_setup"], geometry=cz["source_setup"].split("__")[0],
               pair=cz["source_setup"].split("__")[1].replace("-", "/"), load=case["fe"]["load"],
               k_sigma=cz["k_sigma"], k_G=cz["k_G"], w_mm=cz["w"], lch_over_w=ffm["lch_over_w"],
               sigma_ffm=ffm["sigma_nom_c"], ffm_path=ffm["pid"], ffm_kind=ffm["kind"],
               sigma_onset=sig5, sigma_first_failure=first_fail, sigma_instability=float(s[inst]),
               sigma_separation=sig_sep, sigma_max=float(s[i_max]),
               reserve=(sig_sep / sig5) if np.isfinite(sig_sep) else np.nan,
               separated=bool(len(sep)), completed=r["completed"],
               u_end_over_u_ffm=float(fr.u.iloc[-1] / ffm["delta_onset"]),
               viscous_energy_share=_visc_share(r))
    row["czm_over_ffm"] = row["sigma_onset"] / row["sigma_ffm"]
    curve = pd.DataFrame(dict(name=case["name"], u_over_u_ffm=fr.u / ffm["delta_onset"],
                              sigma=fr.sigma, frac_failed=fr.get("frac_failed", 0.0)))
    return row, curve


def _visc_share(r):
    """(viscous + stabilisation dissipation) / external work at the end of the run."""
    h = r.get("history", {})
    try:
        num = sum(h[k][-1][1] for k in ("ALLVD", "ALLSD") if k in h)
        return float(num / max(h["ALLWK"][-1][1], 1e-30))
    except (KeyError, IndexError, TypeError):
        return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_czm")
    ap.add_argument("--out", default="results_czm")
    ap.add_argument("--fig", default="figures/data")
    a = ap.parse_args()
    rows, curves = [], []
    for c in sorted((ROOT / a.runs).glob("*/case.json")):
        if not (c.parent / "raw" / "czm.json").exists():
            print("   no result:", c.parent.name); continue
        row, cur = analyse(c.parent)
        if row:
            rows.append(row); curves.append(cur)
    if not rows:
        sys.exit("no CZM results")
    S = pd.DataFrame(rows)
    # gain with respect to the flat interface at the same interface parameters and size (sweep)
    flat = S[S.geometry == "flat"].set_index(["pair", "load", "k_sigma", "k_G", "w_mm"])
    def gain(r, col):
        k = (r.pair, r.load, r.k_sigma, r.k_G, r.w_mm)
        return r[col] / flat.loc[k, col] if k in flat.index else np.nan
    S["gain_czm"] = [gain(r, "sigma_onset") for _, r in S.iterrows()]
    S["gain_ffm"] = [gain(r, "sigma_ffm") for _, r in S.iterrows()]
    (ROOT / a.out).mkdir(exist_ok=True)
    S.to_csv(ROOT / a.out / "czm_summary.csv", index=False)
    fig = ROOT / a.fig
    fig.mkdir(parents=True, exist_ok=True)
    pd.concat(curves).to_csv(fig / "fig5a_czm_curves.csv", index=False)
    S.to_csv(fig / "fig5b_czm_summary.csv", index=False)
    ok = S[S.completed | (S.reserve > 1.0)]
    print(f"{len(S)} CZM runs ({int(S.completed.sum())} ran to the end)")
    print(f"CZM onset (5 % non-linearity) / FFM onset: median {S.czm_over_ffm.median():.3f}, "
          f"10-90 %: {S.czm_over_ffm.quantile(0.1):.3f} - {S.czm_over_ffm.quantile(0.9):.3f}")
    print(f"reserve (separation / CZM onset): median {S.reserve.median():.2f}, "
          f"max {S.reserve.max():.2f}; interface fully separated in {int(S.separated.sum())}/{len(S)} runs")
    print(f"viscous dissipation / external work: median {S.viscous_energy_share.median():.2e}")


if __name__ == "__main__":
    main()
