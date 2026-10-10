"""
Finite Fracture Mechanics (coupled stress-energy criterion) kernel.

All inputs are evaluated at a *reference load* (lambda = 1). Linear elasticity
means stresses scale with lambda and energy release rates with lambda**2, so the
critical load factor follows from a single uncracked analysis plus one cracked
analysis per crack length:

    lambda_s(da) = max_{r in [0, da]} 1 / f_s(sigma(r), tau(r))      (stress condition everywhere on the new crack)
    lambda_G(da) = 1 / sqrt(f_G(G_I_inc(da), G_II_inc(da)))            (incremental energy condition)
    lambda_c     = min_da max(lambda_s(da), lambda_G(da))

Units: stresses in MPa, lengths in mm, G in N/mm (= kJ/m^2).
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional, Sequence

import numpy as np

try:
    from scipy.interpolate import PchipInterpolator as _Pchip
except ImportError:          # scipy is optional; falls back to linear log-log interpolation
    _Pchip = None

_BIG = 1e30


@dataclass(frozen=True)
class Criterion:
    """Failure properties of one potential crack plane (interface or bulk)."""
    sigma_c: float                      # tensile strength normal to the crack plane [MPa]
    G_Ic: float                         # mode I toughness [N/mm]
    tau_c: Optional[float] = None       # shear strength [MPa]; None -> normal-stress criterion only
    G_IIc: Optional[float] = None       # mode II toughness [N/mm]; None -> total G compared with G_Ic
    energy_law: str = "linear"          # 'linear' | 'power' | 'bk' | 'total'
    power: float = 1.0                  # exponent for 'power' law
    bk_eta: float = 2.0                 # Benzeggagh-Kenane exponent for 'bk'

    def scaled(self, k_sigma: float = 1.0, k_G: float = 1.0) -> "Criterion":
        """Interface-quality scaling (strengths by k_sigma, toughnesses by k_G)."""
        return replace(
            self,
            sigma_c=self.sigma_c * k_sigma,
            tau_c=None if self.tau_c is None else self.tau_c * k_sigma,
            G_Ic=self.G_Ic * k_G,
            G_IIc=None if self.G_IIc is None else self.G_IIc * k_G,
        )

    # ---- stress condition ------------------------------------------------
    def stress_measure(self, sn, tau):
        """Failure index at unit load; homogeneous of degree 1 in the load."""
        sn = np.asarray(sn, float)
        tau = np.asarray(tau, float)
        sn_pos = np.maximum(sn, 0.0)            # Macaulay bracket: compression does not open a crack
        if self.tau_c is None:
            return sn_pos / self.sigma_c
        return np.sqrt((sn_pos / self.sigma_c) ** 2 + (tau / self.tau_c) ** 2)

    # ---- energy condition ------------------------------------------------
    def energy_index(self, GI, GII):
        """f_G at unit load (degree 2 in the load). lambda_G = 1/sqrt(f_G)."""
        GI = np.maximum(np.asarray(GI, float), 0.0)
        GII = np.maximum(np.asarray(GII, float), 0.0)
        GT = GI + GII
        law = self.energy_law if self.G_IIc is not None else "total"
        if law == "total":
            return GT / self.G_Ic
        if law == "linear":
            return GI / self.G_Ic + GII / self.G_IIc
        if law == "power":
            p = self.power
            return ((GI / self.G_Ic) ** p + (GII / self.G_IIc) ** p) ** (1.0 / p)
        if law == "bk":
            with np.errstate(invalid="ignore", divide="ignore"):
                psi = np.where(GT > 0, GII / GT, 0.0)
            Gc = self.G_Ic + (self.G_IIc - self.G_Ic) * psi ** self.bk_eta
            return GT / Gc
        raise ValueError(f"unknown energy law {law!r}")


@dataclass
class FFMResult:
    lam_c: float            # critical load factor (multiplies the reference load)
    da_c: float             # crack length at onset (per tip for centred cracks)
    lam_s: float            # stress-condition load factor at da_c
    lam_G: float            # energy-condition load factor at da_c
    governing: str          # 'coupled' | 'stress' | 'energy'
    at_lower_bound: bool    # optimum at smallest da -> refine the da grid downwards / mesh
    at_upper_bound: bool    # optimum at largest da  -> extend the path / da grid
    curves: Optional[dict] = None


def _lambda_stress_envelope(branches, grid, crit: Criterion):
    """lambda_s(da): running maximum of the pointwise stress load factor.

    branches: sequence of (r, sn, tau); r = distance from the crack origin (>=0).
    Several branches are used for centred cracks (both sides must satisfy the
    stress condition), their failure indices are combined with max().
    """
    pts = [grid]
    for r, _, _ in branches:
        pts.append(np.asarray(r, float))
    fine = np.unique(np.concatenate(pts))
    fine = fine[fine >= 0.0]
    val = np.zeros_like(fine)
    for r, sn, tau in branches:
        r = np.asarray(r, float)
        order = np.argsort(r)
        v = crit.stress_measure(np.asarray(sn)[order], np.asarray(tau)[order])
        # A point with no stress data on this branch contributes nothing.
        vb = np.interp(fine, r[order], v, left=v[0], right=v[-1])
        val = np.maximum(val, vb)
    # Branch data must cover the crack: points on the new crack are all needed.
    with np.errstate(divide="ignore"):
        lam_pt = np.where(val > 0, 1.0 / val, _BIG)
    # Condition must hold at every point of (0, da]: the load must exceed the
    # largest pointwise lambda. r=0 is excluded only if it is singular (lam->0)
    # which does not change the max, so keep it.
    lam_run = np.maximum.accumulate(lam_pt)
    return np.interp(grid, fine, lam_run)


def _lambda_energy(da, GI, GII, grid, crit: Criterion):
    da = np.asarray(da, float)
    order = np.argsort(da)
    f = crit.energy_index(np.asarray(GI)[order], np.asarray(GII)[order])
    with np.errstate(divide="ignore"):
        lam = np.where(f > 0, 1.0 / np.sqrt(f), _BIG)
    # monotone cubic interpolation in log-log space (G_inc ~ power law near singular points)
    x, y = np.log(da[order]), np.log(lam)
    if _Pchip is not None and len(x) >= 3:
        return np.exp(_Pchip(x, y, extrapolate=True)(np.log(grid)))
    return np.exp(np.interp(np.log(grid), x, y))


def ffm_path(branches: Sequence, da, GI, GII, crit: Criterion,
             n_eval: int = 800, return_curves: bool = False,
             coupled_tol: float = 0.02) -> FFMResult:
    """Critical load factor for one candidate crack path.

    branches : list of (r, sn, tau) stress profiles at reference load, r >= 0
    da       : crack lengths (per tip) for which G_inc was computed
    GI, GII  : incremental energy release rates at reference load [N/mm]
    """
    da = np.asarray(da, float)
    if da.size < 2:
        raise ValueError("need at least two crack lengths")
    grid = np.geomspace(da.min(), da.max(), n_eval)
    lam_s = _lambda_stress_envelope(branches, grid, crit)
    lam_G = _lambda_energy(da, GI, GII, grid, crit)
    lam = np.maximum(lam_s, lam_G)
    i = int(np.argmin(lam))
    at_lo, at_hi = (i == 0), (i == n_eval - 1)
    # local refinement around the coarse optimum
    lo, hi = grid[max(i - 2, 0)], grid[min(i + 2, n_eval - 1)]
    loc = np.geomspace(lo, hi, 400)
    ls_loc = _lambda_stress_envelope(branches, loc, crit)
    lg_loc = _lambda_energy(da, GI, GII, loc, crit)
    l_loc = np.maximum(ls_loc, lg_loc)
    j = int(np.argmin(l_loc))
    lam_c, da_c, ls, lg = l_loc[j], loc[j], ls_loc[j], lg_loc[j]
    if abs(ls - lg) <= coupled_tol * max(ls, lg):
        gov = "coupled"
    else:
        gov = "stress" if ls > lg else "energy"
    curves = dict(da=grid, lam_s=lam_s, lam_G=lam_G) if return_curves else None
    return FFMResult(lam_c=float(lam_c), da_c=float(da_c), lam_s=float(ls),
                     lam_G=float(lg), governing=gov,
                     at_lower_bound=at_lo, at_upper_bound=at_hi, curves=curves)
