"""
Parametric interlocking-interface geometries and candidate FFM crack paths.

Conventions
-----------
* Specimen occupies [0, W] x [-H, H]; material A below the interface, B above.
* The interface polyline runs from (0, 0) to (W, 0); features of A protrude
  upward into B.  Interface arc-length coordinate S increases left -> right.
* Interface frame: t = +S direction, n = left normal (points into B).
* All lengths are in units of the feature width w (FE is run with w = 1 mm; the
  physical size enters only in post-processing, see postprocess.py).

Every candidate crack tip (for every path and every crack length) is written as
a geometric breakpoint, so the uncracked and all cracked models share the same
partitioning and therefore the same mesh; only the seam assignment differs.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, asdict
from typing import List, Optional

import numpy as np

TURN_TOL_DEG = 15.0          # vertices turning more than this are singular corners
SNAP_TOL = 1e-4              # tip points closer than this to a breakpoint are snapped
MIN_FACE_ANGLE_DEG = 15.0    # bulk paths must stay this far from wedge faces


# --------------------------------------------------------------------------
# geometry specification
# --------------------------------------------------------------------------
@dataclass
class GeomSpec:
    shape: str                 # 'flat' | 'rect' | 'U' | 'T' | 'dovetail' | 'sine'
    ar: float = 1.0            # depth / width (d / w)
    alpha_deg: float = 20.0    # dovetail flank angle from vertical
    w: float = 1.0             # neck width (dovetail, T, U, rect); half wavelength (sine)
    pitch: float = 4.0         # cell width in units of w
    n_features: int = 1
    H: float = 8.0             # half height of the specimen in units of w
    head_ratio: float = 2.0    # T: head width / neck width
    neck_frac: float = 0.5     # T: neck depth / total depth
    n_arc: int = 24            # U: segments on the semicircular tip
    n_sine: int = 48           # sine: samples per wavelength
    margin: float = 0.25       # minimum ligament to the cell boundary, units of w

    @property
    def W(self):
        return self.pitch * self.n_features * self.w

    @property
    def d(self):
        return self.ar * self.w

    def tag(self):
        if self.shape == "flat":
            return "flat"
        s = f"{self.shape}_ar{self.ar:g}"
        if self.shape == "dovetail":
            s += f"_a{self.alpha_deg:g}"
        return s


def _feature(g: GeomSpec, xc: float):
    w, d = g.w, g.d
    if g.shape == "rect":
        return [(xc - w / 2, 0), (xc - w / 2, d), (xc + w / 2, d), (xc + w / 2, 0)], w
    if g.shape == "dovetail":
        wt = w + 2 * d * math.tan(math.radians(g.alpha_deg))
        return [(xc - w / 2, 0), (xc - wt / 2, d), (xc + wt / 2, d), (xc + w / 2, 0)], wt
    if g.shape == "T":
        wh, dn = g.head_ratio * w, g.neck_frac * d
        return [(xc - w / 2, 0), (xc - w / 2, dn), (xc - wh / 2, dn), (xc - wh / 2, d),
                (xc + wh / 2, d), (xc + wh / 2, dn), (xc + w / 2, dn), (xc + w / 2, 0)], wh
    if g.shape == "U":
        r = w / 2
        if d < r - 1e-12:
            raise ValueError("U-shape needs ar >= 0.5")
        yc = d - r
        pts = [(xc - w / 2, 0)]
        if yc > 1e-12:
            pts.append((xc - w / 2, yc))
        for th in np.linspace(math.pi, 0, g.n_arc + 1)[1:-1]:
            pts.append((xc + r * math.cos(th), yc + r * math.sin(th)))
        if yc > 1e-12:
            pts.append((xc + w / 2, yc))
        pts.append((xc + w / 2, 0))
        return pts, w
    raise ValueError(g.shape)


def interface_polyline(g: GeomSpec) -> np.ndarray:
    W = g.W
    if g.shape == "flat":
        return np.array([[0.0, 0.0], [W, 0.0]])
    if g.shape == "sine":
        lam = 2 * g.w
        n = max(int(round(W / lam * g.n_sine)), 8)
        x = np.linspace(0, W, n + 1)
        # cosine phase: the interface meets both free edges at a trough (slope 0, 90 deg wedges as
        # in the flat case) and the cell is mirror-symmetric about x = W/2. A sine phase would make
        # the free-edge wedge angles depend on an arbitrary phase choice (artificial singularity).
        y = -0.5 * g.d * np.cos(2 * np.pi * x / lam)
        y[np.abs(y) < 1e-12] = 0.0
        return np.column_stack([x, y])
    pts = [(0.0, 0.0)]
    p = g.pitch * g.w
    for i in range(g.n_features):
        xc = (i + 0.5) * p
        f, width = _feature(g, xc)
        if width / 2 + g.margin * g.w > p / 2 + 1e-12:
            raise ValueError(f"{g.tag()}: feature width {width:.3f} does not fit the cell {p:.3f}")
        pts += f
    pts.append((W, 0.0))
    return np.array(pts, float)


def arc_length(P):
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])


def point_at(P, S, s):
    s = np.clip(s, 0, S[-1])
    i = int(np.clip(np.searchsorted(S, s, side="right") - 1, 0, len(S) - 2))
    t = (s - S[i]) / (S[i + 1] - S[i])
    return P[i] + t * (P[i + 1] - P[i])


def _ang(v):
    return math.degrees(math.atan2(v[1], v[0])) % 360.0


@dataclass
class Singular:
    index: int
    kind: str                 # 'free_edge' | 'corner'
    xy: tuple
    S: float
    wedge_A: float            # interior angle of A at the vertex [deg]
    wedge_B: float
    start_A: float            # CCW start direction of the A wedge [deg]
    start_B: float


def singular_points(P, smooth=False) -> List[Singular]:
    """smooth=True: interior vertices are samples of a smooth curve (sine) -> only free edges."""
    S = arc_length(P)
    out = []
    n = len(P)
    for i in range(n):
        if i == 0 or i == n - 1:
            # free edge: lateral boundary is vertical
            if i == 0:
                phi_next = _ang(P[1] - P[0])
                wB = (90.0 - phi_next) % 360.0
                out.append(Singular(i, "free_edge", tuple(P[i]), S[i], 180.0 - wB, wB, 270.0, phi_next))
            else:
                phi_prev = _ang(P[i - 1] - P[i])
                wB = (phi_prev - 90.0) % 360.0
                out.append(Singular(i, "free_edge", tuple(P[i]), S[i], 180.0 - wB, wB, phi_prev, 90.0))
            continue
        if smooth:
            continue
        phi_prev = _ang(P[i - 1] - P[i])
        phi_next = _ang(P[i + 1] - P[i])
        wB = (phi_prev - phi_next) % 360.0
        wA = 360.0 - wB
        if abs(wB - 180.0) > TURN_TOL_DEG:
            out.append(Singular(i, "corner", tuple(P[i]), S[i], wA, wB, phi_prev, phi_next))
    return out


# --------------------------------------------------------------------------
# candidate paths
# --------------------------------------------------------------------------
@dataclass
class PathSpec:
    pid: str
    kind: str                     # 'interface' | 'bulk'
    mode: str                     # 'edge' | 'centred'
    material: str                 # side whose stresses are used (A/B)
    L: float                      # maximum length (per tip)
    da: List[float] = field(default_factory=list)
    # interface
    S0: Optional[float] = None    # origin (edge) or centre (centred) in interface S
    direction: int = 1            # +1 along S, -1 against S (edge mode)
    s_offset: float = 0.0         # centred: half length of an existing pre-crack
    # bulk
    origin: Optional[tuple] = None
    angle_deg: Optional[float] = None
    note: str = ""
    end: str = ""                 # what limits the path length: feature | obstacle | cap

    def crack_interval(self, da):
        """Interface S-interval(s) of the crack (including pre-crack) for length da."""
        if self.mode == "edge":
            a, b = sorted([self.S0, self.S0 + self.direction * da])
            return a, b
        return self.S0 - self.s_offset - da, self.S0 + self.s_offset + da

    def new_area(self, da):
        return da if self.mode == "edge" else 2.0 * da


def _ray_clip(origin, ang_deg, L, P, box, skip_seg=()):
    """Distance along the ray to the first obstacle (interface or outer box)."""
    o = np.asarray(origin, float)
    e = np.array([math.cos(math.radians(ang_deg)), math.sin(math.radians(ang_deg))])
    best = L
    segs = [(P[i], P[i + 1]) for i in range(len(P) - 1) if i not in skip_seg]
    (x0, y0), (x1, y1) = box
    segs += [(np.array([x0, y0]), np.array([x1, y0])), (np.array([x1, y0]), np.array([x1, y1])),
             (np.array([x1, y1]), np.array([x0, y1])), (np.array([x0, y1]), np.array([x0, y0]))]
    for a, b in segs:
        d = b - a
        M = np.array([[e[0], -d[0]], [e[1], -d[1]]])
        if abs(np.linalg.det(M)) < 1e-14:
            continue
        t, u = np.linalg.solve(M, a - o)
        if t > 1e-9 and -1e-9 <= u <= 1 + 1e-9:
            best = min(best, t)
    return best


def default_da(n=10, da_min=0.005, da_max=0.5):
    return list(np.geomspace(da_min, da_max, n))


def candidate_paths(g: GeomSpec, P, da_grid, L_max=None, bulk_offsets=(-60, -30, 0, 30, 60),
                    bulk_materials=("A", "B"), centred_S=None, precrack=None,
                    extend=None, periodic=False) -> List[PathSpec]:
    """Generate interface and bulk candidate paths from all singular points.

    extend = {"n": .., "min": .., "max": ..}: extended paths (propagation analysis). Interface
    edge paths run to the next singular point (corner / free edge) along the interface,
    capped at extend["max"]; bulk paths to 0.8 x the distance to the next obstacle, capped
    likewise. Each path gets its own geometric da grid from extend["min"] to its length.
    Without extend, the global da_grid and L_max = max(da_grid) are used (onset only)."""
    S = arc_length(P)
    L_max = (extend or {}).get("max") or L_max or max(da_grid)
    # bulk paths may have their own (shorter) cap: a bulk crack longer than ~2w has severed
    # the feature anyway, while interface segments can be longer than 2w (head tops, arcs)
    L_bulk = (extend or {}).get("max_bulk") or L_max
    box = ((0.0, -g.H * g.w), (g.W, g.H * g.w))
    paths: List[PathSpec] = []
    sing = singular_points(P, smooth=(g.shape == "sine"))
    if periodic:
        # interior cell of a periodic array: the cell edges are not free surfaces
        sing = [sp for sp in sing if sp.kind != "free_edge"]
        if not sing and not centred_S and precrack is None:
            centred_S = [(S[-1] / 2, 0.0)]       # flat interface: nucleation in the uniform field
    S_sing = sorted(sp.S for sp in sing)

    def clip_da(L):
        if extend:
            return [float(v) for v in default_da(extend["n"], extend["min"], L)]
        return [d for d in da_grid if d <= L * (1 + 1e-9)]

    def to_next_feature(S0, direction):
        tol = 1e-9 * max(S[-1], 1.0)
        if direction > 0:
            nxt = [v for v in S_sing if v > S0 + tol]
            return (min(nxt) if nxt else S[-1]) - S0
        prv = [v for v in S_sing if v < S0 - tol]
        return S0 - (max(prv) if prv else 0.0)

    d_min = (extend or {}).get("min", min(da_grid))

    # interface paths, edge mode, from every singular point in both directions
    for sp in sing:
        for direction in (+1, -1):
            if precrack is not None:
                continue
            L = min(L_max, (S[-1] - sp.S) if direction > 0 else sp.S)
            end = "feature" if L < L_max else "cap"
            if extend:
                Lf = to_next_feature(sp.S, direction)
                end = "feature" if Lf <= L else end
                L = min(L, Lf)
            if L <= 2 * d_min:
                continue
            paths.append(PathSpec(pid=f"I{sp.index}{'f' if direction > 0 else 'b'}", kind="interface",
                                  mode="edge", material="A", L=L, da=clip_da(L), S0=sp.S,
                                  direction=direction, note=f"{sp.kind} v{sp.index}", end=end))
    # centred interface paths (smooth interfaces, pre-cracked benchmarks)
    for k, (Sc, s_off) in enumerate(centred_S or []):
        Lg = min(Sc - s_off, S[-1] - Sc - s_off)
        L = min(L_max, Lg)
        paths.append(PathSpec(pid=f"C{k}", kind="interface", mode="centred", material="A", L=L,
                              da=clip_da(L), S0=Sc, s_offset=s_off, note="centred",
                              end="feature" if Lg <= L_max else "cap"))
    # bulk paths from re-entrant corners
    for sp in sing:
        if sp.kind != "corner":
            continue
        for mat in bulk_materials:
            wedge = sp.wedge_A if mat == "A" else sp.wedge_B
            start = sp.start_A if mat == "A" else sp.start_B
            if wedge <= 180.0 + 1.0:
                continue
            bis = start + wedge / 2
            for off in bulk_offsets:
                rel = wedge / 2 + off
                if rel < MIN_FACE_ANGLE_DEG or rel > wedge - MIN_FACE_ANGLE_DEG:
                    continue
                ang = (bis + off) % 360.0
                L = 0.8 * _ray_clip(sp.xy, ang, 10 * L_bulk, P, box,
                                    skip_seg=(sp.index - 1, sp.index))
                end = "obstacle" if L < L_bulk else "cap"
                L = min(L, L_bulk)
                if L <= 2 * d_min:
                    continue
                paths.append(PathSpec(pid=f"B{sp.index}{mat}{int(round(ang)):03d}", kind="bulk",
                                      mode="edge", material=mat, L=L, da=clip_da(L),
                                      origin=tuple(sp.xy), angle_deg=ang,
                                      note=f"bulk {mat} from v{sp.index}, wedge {wedge:.0f}", end=end))
    return paths


def mirror_prune(paths: List[PathSpec], P, W, tol=1e-6):
    """Drop paths that are mirror images (x -> W - x) of already kept ones.
    Valid only for loads symmetric about x = W/2; the geometry symmetry is checked here and
    nothing is pruned if the interface polyline is not mirror-symmetric."""
    Pm = np.column_stack([W - P[::-1, 0], P[::-1, 1]])
    if Pm.shape != P.shape or np.max(np.abs(Pm - P)) > 1e-6 * max(W, 1.0):
        return list(paths)
    S = arc_length(P)
    Ltot = S[-1]
    kept, keys = [], []

    def key(p):
        if p.kind == "interface":
            a, b = p.crack_interval(p.L)
            return ("I", p.mode, round(a, 6), round(b, 6))
        ox, oy = p.origin
        return ("B", p.material, round(ox, 6), round(oy, 6), round(p.angle_deg % 360, 3))

    def mkey(p):
        if p.kind == "interface":
            a, b = p.crack_interval(p.L)
            return ("I", p.mode, round(Ltot - b, 6), round(Ltot - a, 6))
        ox, oy = p.origin
        return ("B", p.material, round(W - ox, 6), round(oy, 6), round((180 - p.angle_deg) % 360, 3))

    for p in paths:
        if mkey(p) in keys:
            continue
        keys.append(key(p))
        kept.append(p)
    return kept


# --------------------------------------------------------------------------
# breakpoints and seam selection
# --------------------------------------------------------------------------
def _merge(values, tol):
    values = sorted(values)
    out = []
    for v in values:
        if not out or v - out[-1] > tol:
            out.append(v)
    return out


def _snap(v, ref, tol):
    ref = np.asarray(ref)
    j = int(np.argmin(np.abs(ref - v)))
    return float(ref[j]) if abs(ref[j] - v) <= tol else float(v)


def build_case(g: GeomSpec, paths: List[PathSpec], P, precrack=None):
    """Assemble breakpoints, seam edge points and path metadata for the FE driver."""
    S = arc_length(P)
    tol = SNAP_TOL * g.w
    # ---- interface breakpoints
    bp = list(S)
    if precrack is not None:
        bp += list(precrack)
        # mirror the crack-length breakpoints behind each pre-crack tip: the tip of the
        # uncracked (pre-cracked) model must be meshed as finely as the tips of the cracked
        # models, otherwise the difference of the two discretisation errors appears as a
        # spurious constant offset in Delta(ALLSE) that dominates the smallest increments
        half = 0.5 * (precrack[1] - precrack[0])
        for d in sorted({float(x) for p in paths for x in p.da}):
            if d < 0.9 * half:
                bp += [precrack[0] + d, precrack[1] - d]
    for p in paths:
        if p.kind != "interface":
            continue
        for da in p.da:
            bp += list(p.crack_interval(da))
    bp = [b for b in bp if -1e-12 <= b <= S[-1] + 1e-12]
    bp = _merge(bp, tol)
    bp_arr = np.array(bp)

    def seam_points(a, b):
        a, b = _snap(a, bp_arr, tol), _snap(b, bp_arr, tol)
        inside = bp_arr[(bp_arr >= a - 1e-12) & (bp_arr <= b + 1e-12)]
        mids = 0.5 * (inside[:-1] + inside[1:])
        return [point_at(P, S, m).tolist() for m in mids], (a, b)

    pre = None
    if precrack is not None:
        pts, (a, b) = seam_points(*precrack)
        pre = dict(S0=a, S1=b, seam_points=pts)

    # crossings between straight bulk lines become partition vertices in Abaqus; they must
    # also be breakpoints of both lines, otherwise a crack running through a crossing is
    # only partly seamed (one of the two sub-edges is missed)
    bulk = [p for p in paths if p.kind == "bulk"]
    extra = {p.pid: [] for p in bulk}
    for i in range(len(bulk)):
        for j in range(i + 1, len(bulk)):
            pi, pj = bulk[i], bulk[j]
            oi, oj = np.array(pi.origin), np.array(pj.origin)
            ei = np.array([math.cos(math.radians(pi.angle_deg)), math.sin(math.radians(pi.angle_deg))])
            ej = np.array([math.cos(math.radians(pj.angle_deg)), math.sin(math.radians(pj.angle_deg))])
            den = ei[0] * (-ej[1]) - ei[1] * (-ej[0])
            if abs(den) < 1e-9:
                continue                      # parallel
            r = oj - oi
            si = (r[0] * (-ej[1]) - r[1] * (-ej[0])) / den
            sj = (ei[0] * r[1] - ei[1] * r[0]) / den
            if tol < si < pi.L - tol and tol < sj < pj.L - tol:
                extra[pi.pid].append(float(si))
                extra[pj.pid].append(float(sj))

    bulk_lines, out_paths = [], []
    for p in paths:
        d = asdict(p)
        cracks = []
        if p.kind == "interface":
            for da in p.da:
                a, b = p.crack_interval(da)
                pts, (a2, b2) = seam_points(a, b)
                cracks.append(dict(da=float(da), S_interval=[a2, b2], seam_points=pts))
        else:
            o = np.array(p.origin)
            e = np.array([math.cos(math.radians(p.angle_deg)), math.sin(math.radians(p.angle_deg))])
            ss = _merge([0.0] + list(p.da) + [p.L] + extra.get(p.pid, []), tol)
            line_pts = [(o + s * e).tolist() for s in ss]
            bulk_lines.append(dict(id=p.pid, material=p.material, points=line_pts,
                                   origin=list(o), angle_deg=p.angle_deg, L=p.L))
            for da in p.da:
                inside = [s for s in ss if s <= da + 1e-12]
                mids = [0.5 * (inside[k] + inside[k + 1]) for k in range(len(inside) - 1)]
                cracks.append(dict(da=float(da), seam_points=[(o + m * e).tolist() for m in mids],
                                   line=[list(o), (o + da * e).tolist()]))
        d["cracks"] = cracks
        out_paths.append(d)

    # points guaranteed to lie inside A and B (far from the interface)
    H = g.H * g.w
    return dict(
        geometry=asdict(g), tag=g.tag(), W=g.W, H=H,
        interface_points=P.tolist(),
        interface_breakpoints=[point_at(P, S, b).tolist() + [b] for b in bp],
        bulk_lines=bulk_lines,
        point_in_A=[g.W / 2, -0.5 * H], point_in_B=[g.W / 2, 0.75 * H],
        singular=[asdict(s) for s in singular_points(P, smooth=(g.shape == "sine"))],
        precrack=pre,
        paths=out_paths,
    )


def sine_centres(g, P):
    """Centres of centred interface paths on a smooth (sine) interface: every quarter
    wavelength (crests, troughs and inflection points), interior only. Returns [(S, 0.0)]."""
    S = arc_length(P)
    lam = 2 * g.w
    out = []
    for k in range(1, int(round(4 * g.W / lam))):
        x = k * lam / 4
        i = int(np.argmin(np.abs(P[:, 0] - x)))
        out.append((float(S[i]), 0.0))
    return out
