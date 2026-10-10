"""Unit tests (run: python -m pytest -q tests)."""
import sys, pathlib, math, json
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import pytest
from ffm.abq_utils import Polyline, nodes_on_polyline, node_to_elements, seam_jumps, side_by_material, \
    side_by_centroid, sum_nodal_forces, left_of_line_fn
from ffm.kernel import Criterion, ffm_path
from ffm.materials import dundurs
from ffm.geometry import GeomSpec, interface_polyline, candidate_paths, mirror_prune, build_case, default_da, singular_points


def test_polyline_projection():
    pl = Polyline([(0, 0), (1, 0), (1, 1)])
    S, d, seg = pl.project(1.0, 0.5)
    assert abs(S - 1.5) < 1e-12 and d < 1e-12 and seg == 1
    S, d, seg = pl.project(0.3, 0.2)
    assert abs(S - 0.3) < 1e-12 and abs(d - 0.2) < 1e-12


def test_seam_pairing_interface():
    # two quads sharing an interface y=0 with duplicated nodes 10/11 (A side) and 20/21 (B side)
    coords = {10: (0.0, 0.0), 11: (1.0, 0.0), 20: (0.0, 0.0), 21: (1.0, 0.0), 1: (0, -1), 2: (1, -1), 3: (0, 1), 4: (1, 1)}
    conn = {100: (1, 2, 11, 10), 200: (20, 21, 4, 3)}
    u = {10: (0, 0), 11: (0, 0), 20: (0.01, 0.1), 21: (0.02, 0.2), 1: (0, 0), 2: (0, 0), 3: (0, 0), 4: (0, 0)}
    pl = Polyline([(0, 0), (1, 0)])
    on = nodes_on_polyline(coords, pl, 1e-9)
    n2e = node_to_elements(conn, [r[0] for r in on])
    J = seam_jumps(on, coords, u, n2e, side_by_material(n2e, {200}), 1e-9)
    assert len(J) == 2
    assert abs(J[0][3] - 0.1) < 1e-12 and abs(J[1][3] - 0.2) < 1e-12   # u(B) - u(A)
    side = side_by_centroid(pl, coords, conn, n2e)
    assert side(20, 0.0, 0) == 1 and side(10, 0.0, 0) == -1


def test_nodal_force_side_sum():
    coords = {1: (0, 0), 2: (1, 0), 3: (0, 1), 4: (0, -1)}
    conn = {7: (1, 2, 3), 8: (1, 4, 2)}
    vals = [(7, 1, (0.0, 2.0)), (8, 1, (0.0, -2.0)), (7, 2, (0.0, 1.0))]
    ok = left_of_line_fn((0, 0), (1, 0), coords, conn)
    f = sum_nodal_forces(vals, {1, 2}, ok)
    assert f[1] == [0.0, 2.0] and f[2] == [0.0, 1.0]


def test_mixed_mode_laws():
    c_lin = Criterion(sigma_c=30, tau_c=20, G_Ic=1.0, G_IIc=2.0, energy_law="linear")
    c_pow = Criterion(sigma_c=30, tau_c=20, G_Ic=1.0, G_IIc=2.0, energy_law="power", power=1.0)
    assert np.isclose(c_lin.energy_index(0.3, 0.4), c_pow.energy_index(0.3, 0.4))
    c_bk = Criterion(sigma_c=30, tau_c=20, G_Ic=1.0, G_IIc=2.0, energy_law="bk", bk_eta=2.0)
    assert np.isclose(c_bk.energy_index(1.0, 0.0), 1.0) and np.isclose(c_bk.energy_index(0.0, 2.0), 1.0)
    # compression does not count in the stress criterion
    assert c_lin.stress_measure(-50.0, 0.0) == 0.0


def test_scaling_homogeneity():
    """lambda_c must scale as 1/load: doubling the reference stresses and 4x G halves lambda."""
    crit = Criterion(sigma_c=40, G_Ic=1.0)
    r = np.linspace(0, 2, 200); sn = 30 / np.sqrt(r + 0.01)
    da = np.geomspace(0.01, 1, 10); G = 0.2 * da ** 0.3
    a = ffm_path([(r, sn, 0 * r)], da, G, 0 * G, crit).lam_c
    b = ffm_path([(r, 2 * sn, 0 * r)], da, 4 * G, 0 * G, crit).lam_c
    assert np.isclose(a, 2 * b, rtol=1e-6)


def test_dundurs_identical():
    al, be, eps = dundurs(3000, 0.35, 3000, 0.35)
    assert abs(al) < 1e-12 and abs(be) < 1e-12 and abs(eps) < 1e-12


@pytest.mark.parametrize("shape,ar,alpha", [("flat", 1, 20), ("U", 0.5, 20), ("T", 2, 20), ("dovetail", 2, 30), ("sine", 2, 20)])
def test_geometry_cases(shape, ar, alpha):
    g = GeomSpec(shape=shape, ar=ar, alpha_deg=alpha)
    P = interface_polyline(g)
    assert np.all(np.diff(P[:, 0]) > -1e-12) or shape in ("T", "dovetail", "U")
    paths = mirror_prune(candidate_paths(g, P, default_da()), P, g.W)
    case = build_case(g, paths, P)
    # every crack has at least one seam edge point, lying on the interface or on its bulk line
    pl = Polyline(case["interface_points"])
    for p in case["paths"]:
        for c in p["cracks"]:
            assert c["seam_points"]
            if p["kind"] == "interface":
                for x, y in c["seam_points"]:
                    assert pl.project(x, y)[1] < 1e-9
    # the points used to pick the A/B faces are on the correct side
    assert case["point_in_A"][1] < P[:, 1].min() and case["point_in_B"][1] > P[:, 1].max()


def test_dovetail_wedges():
    g = GeomSpec(shape="dovetail", ar=1, alpha_deg=20)
    sp = singular_points(interface_polyline(g))
    corners = [s for s in sp if s.kind == "corner"]
    assert np.isclose(corners[0].wedge_A, 290) and np.isclose(corners[1].wedge_B, 290)


def test_sine_symmetric_and_pruning_guard():
    from ffm.geometry import GeomSpec, interface_polyline, candidate_paths, mirror_prune, default_da, sine_centres
    g = GeomSpec("sine", ar=1.0)
    P = interface_polyline(g)
    # cosine phase: horizontal tangent at both free edges, mirror symmetric cell
    assert abs(P[1, 1] - P[0, 1]) < 0.02 * g.d and np.allclose(P[:, 1], P[::-1, 1])
    paths = candidate_paths(g, P, default_da(), centred_S=sine_centres(g, P))
    assert len(mirror_prune(paths, P, g.W)) < len(paths)
    # an asymmetric polyline must not be pruned
    Q = P.copy(); Q[:, 1] = 0.5 * np.sin(np.pi * Q[:, 0])
    assert len(mirror_prune(paths, Q, g.W)) == len(paths)


def test_closure_fallback_when_energy_unresolved():
    from ffm.postprocess import PathData
    da = np.array([0.01, 0.02, 0.04])
    pd = PathData(pid="x", kind="interface", mode="edge", material="A", branches=[], da=da,
                  G_energy=np.array([0.0, 2.0, 3.0]), GI_cc=np.array([0.9, 1.5, 2.4]),
                  GII_cc=np.array([0.1, 0.5, 0.6]), min_opening=np.zeros(3),
                  energy_ok=np.array([False, True, True]))
    GI, GII = pd.modes()
    assert np.allclose(GI + GII, [1.0, 2.0, 3.0])        # closure used where Delta ALLSE = 0
    assert np.isnan(pd.cc_ratio[0]) and np.isclose(pd.cc_ratio[1], 1.0)
