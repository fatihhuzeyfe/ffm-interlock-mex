"""Elastic helpers and dimensionless groups."""
from __future__ import annotations

import json
import math
from pathlib import Path


def E_prime(E, nu, plane="strain"):
    return E / (1.0 - nu ** 2) if plane == "strain" else E


def kolosov(nu, plane="strain"):
    return 3.0 - 4.0 * nu if plane == "strain" else (3.0 - nu) / (1.0 + nu)


def dundurs(E1, nu1, E2, nu2, plane="strain"):
    """Dundurs alpha, beta and the oscillation index eps for material 1 / material 2."""
    mu1, mu2 = E1 / (2 * (1 + nu1)), E2 / (2 * (1 + nu2))
    k1, k2 = kolosov(nu1, plane), kolosov(nu2, plane)
    G = mu1 / mu2
    den = G * (k2 + 1) + (k1 + 1)
    alpha = (G * (k2 + 1) - (k1 + 1)) / den
    beta = (G * (k2 - 1) - (k1 - 1)) / den
    eps = math.log((1 - beta) / (1 + beta)) / (2 * math.pi)
    return alpha, beta, eps


def E_star(E1, nu1, E2, nu2, plane="strain"):
    """Effective bimaterial modulus 2/E* = 1/E1' + 1/E2' (Hutchinson & Suo)."""
    return 2.0 / (1.0 / E_prime(E1, nu1, plane) + 1.0 / E_prime(E2, nu2, plane))


def l_ch(E, G_c, sigma_c):
    """Irwin-type characteristic length E G_c / sigma_c^2 [mm]."""
    return E * G_c / sigma_c ** 2


def brittleness(size, E, G_c, sigma_c):
    """mu = sigma_c^2 * l / (E * G_c) = l / l_ch."""
    return size / l_ch(E, G_c, sigma_c)


def load_materials(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    placeholders = [k for k, v in data.get("materials", {}).items()
                    if isinstance(v, dict) and v.get("_status") == "PLACEHOLDER"]
    placeholders += [k for k, v in data.get("interfaces", {}).items()
                     if isinstance(v, dict) and v.get("_status") == "PLACEHOLDER"]
    if placeholders:
        print("WARNING: placeholder properties still in use for:", ", ".join(placeholders))
    return data


# ----------------------------------------------------------------------------- orthotropy
# Material axes (MEX, transversely isotropic about the raster axis):
#   1 = raster / road direction, 2 = in-layer transverse (road-to-road, intralayer bond),
#   3 = build direction Z (interlayer bond). 2D model plane = layer plane (1-2), plane strain in 3.
# Strengths: X along 1, Y transverse (= Z, transverse isotropy), S in-plane shear (1-2).
# Toughness of a bulk crack: G_L for a crack running ALONG the roads (separating them),
# G_T for a crack running ACROSS the roads.

def tsai_hill_uniaxial(beta_deg, X, Y, S):
    """Uniaxial tensile strength in a direction at angle beta to axis 1 (Tsai-Hill)."""
    b = math.radians(beta_deg)
    c2, s2 = math.cos(b) ** 2, math.sin(b) ** 2
    inv = c2 * c2 / X ** 2 + s2 * s2 / Y ** 2 + (1.0 / S ** 2 - 1.0 / X ** 2) * c2 * s2
    return 1.0 / math.sqrt(inv)


def toughness_orient(psi_deg, G_L, G_T):
    """Bulk toughness for a crack at angle psi to the roads (0 = along the roads)."""
    p = math.radians(psi_deg)
    return G_L * math.cos(p) ** 2 + G_T * math.sin(p) ** 2


def bulk_props_for_crack(mat, crack_angle_deg, raster_deg):
    """sigma_c (normal to the crack plane) and G_c for a straight bulk crack at
    crack_angle_deg (global) in a material whose roads run at raster_deg (global)."""
    psi = (crack_angle_deg - raster_deg) % 180.0
    beta = psi + 90.0                      # crack-plane normal relative to the roads
    if "X" in mat:
        sc = tsai_hill_uniaxial(beta, mat["X"], mat["Y"], mat["S"])
        gc = toughness_orient(psi, mat["G_L"], mat["G_T"])
    else:                                  # isotropic entry
        sc, gc = mat["sigma_c"], mat["G_c"]
    return sc, gc


def interface_props(mats, A, B):
    """Interface entry for pair A/B. Explicit values win; with "rule": "weaker_intralayer"
    the reference interface (k = 1) is the intralayer (road-to-road) bond of the weaker
    material: sigma_c = min(Y), tau_c = min(S), G_Ic = min(G_L), G_IIc = ratio * G_Ic."""
    I = mats["interfaces"].get(f"{A}/{B}") or mats["interfaces"].get(f"{B}/{A}")
    if I is None:
        raise KeyError(f"no interface entry '{A}/{B}'")
    out = dict(I)
    if I.get("rule") == "weaker_intralayer":
        ma, mb = mats["materials"][A], mats["materials"][B]
        out.setdefault("sigma_c", min(ma["Y"], mb["Y"]))
        out.setdefault("tau_c", min(ma["S"], mb["S"]))
        out.setdefault("G_Ic", min(ma["G_L"], mb["G_L"]))
        out.setdefault("G_IIc", I.get("GIIc_over_GIc", 2.0) * out["G_Ic"])
    return out


def ortho_fe_entry(mat):
    """Elastic block written to case.json for the Abaqus driver."""
    keys = ("E1", "E2", "E3", "nu12", "nu13", "nu23", "G12", "G13", "G23")
    if all(k in mat for k in keys):
        return {k: float(mat[k]) for k in keys}
    return None


def ortho_check(mat):
    """Positive definiteness of the orthotropic compliance (raises ValueError)."""
    E1, E2, E3 = mat["E1"], mat["E2"], mat["E3"]
    n12, n13, n23 = mat["nu12"], mat["nu13"], mat["nu23"]
    n21, n31, n32 = n12 * E2 / E1, n13 * E3 / E1, n23 * E3 / E2
    det = 1 - n12 * n21 - n23 * n32 - n13 * n31 - 2 * n21 * n32 * n13
    if min(E1, E2, E3, mat["G12"], mat["G13"], mat["G23"]) <= 0 or det <= 0:
        raise ValueError(f"orthotropic constants not positive definite (det = {det:.3f})")
    return det
