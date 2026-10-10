"""
Raw FE output (from abaqus/build_and_run.py or any other solver writing the same
JSON layout) -> per-path stress profiles and incremental energy release rates ->
FFM critical loads.

Energy release rates are taken in magnitude from the strain-energy difference
(U0 - U)/A_new (displacement control) and split into modes with the finite
crack-closure integral

    G_I  = 1/(2 A_new) * int sigma_nn^0 [u_n] ds,   G_II = 1/(2 A_new) * int tau^0 [u_t] ds

(sigma^0: uncracked traction on the future crack line, [u]: jump in the cracked
state). The ratio closure/energy is reported as a mesh-consistency check.

Size scaling: FE is run with feature width w = 1 mm. For a physical width
w_phys = kappa * 1 mm, stresses at homologous points are unchanged while lengths
and G scale with kappa.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .geometry import arc_length
from .kernel import Criterion, ffm_path


def _frame_from_polyline(P, seg):
    d = P[seg + 1] - P[seg]
    t = d / np.linalg.norm(d)
    return t, np.array([-t[1], t[0]])


def traction(sxx, syy, sxy, t, n):
    sn = sxx * n[0] ** 2 + syy * n[1] ** 2 + 2 * sxy * n[0] * n[1]
    tau = sxx * t[0] * n[0] + syy * t[1] * n[1] + sxy * (t[0] * n[1] + t[1] * n[0])
    return sn, tau


# ALLSE is stored in single precision (~6e-8 relative): energy differences smaller than this
# fraction of U0 are not trusted and G is taken from the crack-closure integral instead
ENERGY_REL_TOL = 1e-5


@dataclass
class PathData:
    pid: str
    kind: str
    mode: str
    material: str
    branches: list                 # [(r, sn, tau), ...] at reference load
    da: np.ndarray
    G_energy: np.ndarray           # from Delta ALLSE
    GI_cc: np.ndarray              # crack closure
    GII_cc: np.ndarray
    min_opening: np.ndarray        # most negative [u_n] on the crack (interpenetration check)
    meta: dict = field(default_factory=dict)
    energy_ok: np.ndarray = None   # Delta ALLSE resolvable in the (single precision) ODB
    pen_I: np.ndarray = None       # part of G_I coming from overlapping (penetrating) nodes

    def __post_init__(self):
        if self.energy_ok is None:
            self.energy_ok = np.ones_like(self.da, dtype=bool)
        if self.pen_I is None:
            self.pen_I = np.zeros_like(self.da, dtype=float)

    @property
    def cc_ratio(self):
        with np.errstate(invalid="ignore", divide="ignore"):
            r = (self.GI_cc + self.GII_cc) / self.G_energy
        return np.where(self.energy_ok, r, np.nan)

    def modes(self):
        """G_I, G_II: magnitude from energy (closure where Delta ALLSE is below the output
        precision), partition from closure."""
        tot = self.GI_cc + self.GII_cc
        with np.errstate(invalid="ignore", divide="ignore"):
            phi = np.where(tot > 0, self.GII_cc / tot, 0.0)
        G = np.where(self.energy_ok, self.G_energy, np.maximum(tot, 0.0))
        return G * (1 - phi), G * phi


def crack_extent_check(path, crack, raw, n_per_edge=4):
    """Do the crack-face nodes recorded in a cracked result span exactly this crack?
    Returns (ok, kind, message); kind is 'outside' (nodes beyond the crack: result belongs to a
    different geometry), 'short' (crack faces end before the defined tip, e.g. an unseamed
    sub-edge) or ''. Allowed gap at a tip: one element, (da_k - da_k-1) / n_per_edge."""
    pos = [j[0] for j in raw.get("jumps", [])]
    if not pos:
        return False, "short", "no crack-face nodes"
    da = crack["da"]
    das = sorted(c["da"] for c in path["cracks"])
    k = min(range(len(das)), key=lambda i: abs(das[i] - da))
    seg = das[k] - (das[k - 1] if k > 0 else 0.0)
    gap_ok = 1.5 * seg / n_per_edge + 0.02 * da + 1e-6
    out_tol = 0.02 * da + 1e-5
    if path["kind"] == "bulk":
        if max(pos) > da + out_tol or min(pos) < -out_tol:
            return False, "outside", f"beyond crack (max {max(pos):.4f} > da {da:.4f})"
        if da - max(pos) > gap_ok:
            return False, "short", f"crack shorter than defined ({max(pos):.4f} vs {da:.4f})"
        return True, "", ""
    lo, hi = crack["S_interval"]
    if max(pos) > hi + out_tol or min(pos) < lo - out_tol:
        return False, "outside", f"outside crack interval [{lo:.4f}, {hi:.4f}]"
    if path["mode"] == "centred":
        ends = [min(pos) - lo, hi - max(pos)]
    elif path.get("direction", 1) > 0:
        ends = [hi - max(pos)]
    else:
        ends = [min(pos) - lo]
    if max(ends) > gap_ok:
        return False, "short", f"crack shorter than defined (gap {max(ends):.4f}, da {da:.4f})"
    return True, "", ""


def _drop_origin(r, sn, tau, tol=1e-9):
    r, sn, tau = np.asarray(r), np.asarray(sn), np.asarray(tau)
    keep = r > tol
    if keep.sum() < 2:
        return r, sn, tau
    return r[keep], sn[keep], tau[keep]


class Setup:
    """One FE setup directory: case.json + raw/*.json."""

    def __init__(self, directory, raw_subdir="raw"):
        self.dir = Path(directory)
        self.raw = self.dir / raw_subdir
        self.case = json.loads((self.dir / "case.json").read_text())
        self.P = np.array(self.case["interface_points"], float)
        self.S = arc_length(self.P)
        self.unc = json.loads((self.raw / "uncracked.json").read_text())
        self._prep_interface()

    # -------------------------------------------------------------- uncracked
    def _prep_interface(self):
        rows = self.unc["interface"]
        S, sn, tau, SF, F = [], [], [], [], []
        for row in rows:
            s_, seg, x, y, sA, sB = row[:6]
            SF.append(s_); F.append(row[6] if len(row) > 6 else (0.0, 0.0))
            src = sA if sA[0] is not None else sB
            if src[0] is None:
                continue
            t, n = _frame_from_polyline(self.P, int(seg))
            a, b = traction(*src, t, n)
            S.append(s_); sn.append(a); tau.append(b)
        o = np.argsort(S)
        self.if_S = np.array(S)[o]
        self.if_sn = np.array(sn)[o]
        self.if_tau = np.array(tau)[o]
        o = np.argsort(SF)
        self.if_FS = np.array(SF)[o]
        self.if_F = np.array(F, float).reshape(-1, 2)[o]

    def _bulk_profile(self, pid):
        rows = np.array(self.unc["bulk"][pid], float)
        bl = next(b for b in self.case["bulk_lines"] if b["id"] == pid)
        ang = np.radians(bl["angle_deg"])
        t = np.array([np.cos(ang), np.sin(ang)]); n = np.array([-t[1], t[0]])
        s = rows[:, 0]
        sn, tau = traction(rows[:, 2], rows[:, 3], rows[:, 4], t, n)
        o = np.argsort(s)
        F = rows[:, 5:7] if rows.shape[1] >= 7 else np.zeros((len(s), 2))
        return s[o], sn[o], tau[o], t, n, F[o]

    def branches(self, path):
        """(r, sigma_n, tau) along the candidate crack, r measured from the crack origin.
        The origin node itself (r = 0) is excluded: it is the singular point (corner, crack
        tip) where averaged nodal stresses are meaningless - e.g. at a convex corner of one
        material they come out far too LOW and would block the stress criterion for the
        whole path (running maximum). The stress condition is required on (0, da]."""
        return [_drop_origin(*b) for b in self._branches_raw(path)]

    def _branches_raw(self, path):
        if path["kind"] == "bulk":
            s, sn, tau, _, _, _ = self._bulk_profile(path["pid"])
            return [(s, sn, tau)]
        S0 = path["S0"]
        if path["mode"] == "edge":
            if path["direction"] > 0:
                m = self.if_S >= S0 - 1e-12
                return [(self.if_S[m] - S0, self.if_sn[m], self.if_tau[m])]
            m = self.if_S <= S0 + 1e-12
            return [(S0 - self.if_S[m], self.if_sn[m], self.if_tau[m])]
        off = path["s_offset"]
        mf = self.if_S >= S0 + off - 1e-12
        mb = self.if_S <= S0 - off + 1e-12
        return [(self.if_S[mf] - (S0 + off), self.if_sn[mf], self.if_tau[mf]),
                ((S0 - off) - self.if_S[mb], self.if_sn[mb], self.if_tau[mb])]

    # -------------------------------------------------------------- cracked
    def _closure(self, path, crack, raw):
        """Discrete finite crack closure: W = 1/2 sum_i F_i^0 . [u]_i over the nodes released
        by the new crack. F^0: nodal forces transmitted across the line in the uncracked
        model (NFORC summed over the + side), [u]: jumps in the cracked model. With identical
        meshes this equals U0 - U exactly (up to the sign convention of NFORC, which is
        calibrated per path in path_data). Returns (W_I, W_II, min_opening_ratio)."""
        J = np.array(raw["jumps"], float) if raw["jumps"] else np.zeros((0, 4))
        if path["kind"] == "interface":
            lo, hi = crack["S_interval"]
            Fs, F = self.if_FS, self.if_F
            frames = [_frame_from_polyline(self.P, int(seg)) for seg in J[:, 1]]
        else:
            s0, _, _, t, n, F = self._bulk_profile(path["pid"])
            Fs = s0
            lo, hi = 0.0, crack["da"]
            frames = [(t, n)] * len(J)
        tol = 1e-7 * max(1.0, float(np.max(np.abs(Fs))) if len(Fs) else 1.0)
        WI = WII = Wpen = 0.0
        un_all, du_all = [], []
        for (pos, seg, dux, duy), (t, n) in zip(J, frames):
            if pos < lo - tol or pos > hi + tol:
                continue
            j = int(np.argmin(np.abs(Fs - pos)))
            if abs(Fs[j] - pos) > 1e3 * tol:
                continue                      # no matching uncracked node (should not happen)
            f = F[j]
            du = np.array([dux, duy])
            WI += 0.5 * (f @ n) * (du @ n)
            WII += 0.5 * (f @ t) * (du @ t)
            if du @ n < 0:                     # faces overlap here (no contact in the model)
                Wpen += 0.5 * (f @ n) * (du @ n)
            un_all.append(du @ n)
            du_all.append(float(np.hypot(*du)))
        if un_all:
            # overlap relative to the largest TOTAL face displacement (a mode II crack has tiny
            # normal openings, so normalising by them alone would flag negligible overlaps)
            scale = max(du_all) or 1.0
            mn = float(min(un_all)) / scale
        else:
            mn = np.nan
        return WI, WII, mn, Wpen

    def path_data(self, path) -> Optional[PathData]:
        U0 = self.unc["ALLSE"]
        da, Ge, GI, GII, mo, ok, pen = [], [], [], [], [], [], []
        n_excl = 0
        npe = self.case["fe"]["mesh"].get("n_per_edge", 4)
        for k, crack in enumerate(path["cracks"]):
            f = self.raw / f"{path['pid']}__{k:02d}.json"
            if not f.exists():
                continue
            raw = json.loads(f.read_text())
            if not crack_extent_check(path, crack, raw, npe)[0]:
                n_excl += 1                    # FE crack differs from the definition: not used
                continue
            A = crack["da"] if path["mode"] == "edge" else 2 * crack["da"]
            WI, WII, mn, Wpen = self._closure(path, crack, raw)
            dU = U0 - raw["ALLSE"]
            da.append(crack["da"]); Ge.append(dU / A); ok.append(abs(dU) > ENERGY_REL_TOL * abs(U0))
            GI.append(WI / A); GII.append(WII / A); mo.append(mn); pen.append(Wpen / A)
        if len(da) < 2:
            return None
        GI, GII, Ge, ok = np.array(GI), np.array(GII), np.array(Ge), np.array(ok)
        # NFORC sign convention: choose the sign that makes closure agree with Delta ALLSE
        # (resolvable increments only); without any, the sign that makes closure positive
        tot = GI + GII
        if ok.any():
            sgn = 1.0 if np.nanmedian(tot[ok] / Ge[ok]) >= 0 else -1.0
        else:
            sgn = 1.0 if np.nansum(tot) >= 0 else -1.0
        return PathData(pid=path["pid"], kind=path["kind"], mode=path["mode"],
                        material=path["material"], branches=self.branches(path),
                        da=np.array(da), G_energy=Ge, GI_cc=sgn * GI, GII_cc=sgn * GII,
                        min_opening=np.array(mo), energy_ok=ok, pen_I=sgn * np.array(pen),
                        meta=dict(note=path.get("note", ""), nforc_sign=sgn, angle_deg=path.get("angle_deg"),
                                  n_excluded=n_excl,
                                  end=path.get("end", ""),
                                  n_closure_only=int((~ok).sum())))

    def all_paths(self) -> List[PathData]:
        out = []
        for p in self.case["paths"]:
            d = self.path_data(p)
            if d is not None:
                out.append(d)
        return out


def evaluate(setup: Setup, crit_interface: Criterion, crit_bulk: Dict[str, Criterion],
             kappa: float = 1.0, paths: Optional[List[PathData]] = None):
    """FFM on every path of a setup. Returns list of dicts (one per path), sorted by lambda_c."""
    paths = paths if paths is not None else setup.all_paths()
    W = setup.case["W"]
    F_ref = setup.unc["F"]
    rows = []
    for pd in paths:
        if pd.kind == "interface":
            crit = crit_interface
        else:
            crit = crit_bulk[pd.material]
            if callable(crit):                 # orientation-dependent bulk criterion
                crit = crit(pd.meta.get("angle_deg") or 0.0)
        GI, GII = pd.modes()
        br = [(r * kappa, sn, tau) for r, sn, tau in pd.branches]
        res = ffm_path(br, pd.da * kappa, GI * kappa, GII * kappa, crit)
        # overlap diagnostics at the onset length; first-order estimate without the spurious
        # mode I energy of overlapping nodes (a frictionless-contact proxy)
        x = res.da_c / kappa
        tot = np.maximum(pd.GI_cc + pd.GII_cc, 1e-300)
        pen_share = float(np.interp(x, pd.da, pd.pen_I / tot))
        overlap = float(np.interp(x, pd.da, np.nan_to_num(pd.min_opening, nan=0.0)))
        lam_np = res.lam_c
        if pen_share > 0.01:
            GIc_ = np.maximum(GI - pd.pen_I * GI / np.maximum(pd.GI_cc, 1e-300), 0.0)
            lam_np = ffm_path(br, pd.da * kappa, GIc_ * kappa, GII * kappa, crit).lam_c
        rows.append(dict(
            setup=setup.case["name"], pid=pd.pid, kind=pd.kind, material=pd.material, note=pd.meta["note"],
            lam_c=res.lam_c, sigma_nom_c=res.lam_c * F_ref / W, da_c=res.da_c,
            governing=res.governing, at_lower=res.at_lower_bound,
            # optimum beyond the path: if the path ends at a geometric feature this is physical
            # (the whole segment debonds at onset -> sequential analysis), otherwise extend da
            at_upper=res.at_upper_bound and pd.meta.get("end") not in ("feature", "obstacle"),
            segment_exhausted=res.at_upper_bound and pd.meta.get("end") in ("feature", "obstacle"),
            cc_ratio_median=(float(np.nanmedian(pd.cc_ratio)) if np.isfinite(pd.cc_ratio).any() else float('nan')),
            interpenetration=bool(overlap < -0.01 and pen_share > 0.01),   # -> contact matters
            overlap_ratio=overlap, pen_share=pen_share, lam_c_nopen=lam_np,
            n_closure_only=pd.meta.get("n_closure_only", 0),
            mode_mix_at_c=float(np.interp(res.da_c / kappa, pd.da, GII / np.maximum(GI + GII, 1e-300))),
        ))
    rows.sort(key=lambda r: r["lam_c"])
    return rows
