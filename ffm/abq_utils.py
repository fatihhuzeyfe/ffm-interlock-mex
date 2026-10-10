# -*- coding: utf-8 -*-
"""
Pure-Python geometry/extraction helpers that run inside the Abaqus Python
interpreter (2.7 or 3.x). No numpy, no f-strings, no type hints on purpose.
They are unit-tested outside Abaqus (tests/test_abq_utils.py).
"""
from __future__ import division, print_function

import math


class Polyline(object):
    def __init__(self, pts):
        self.pts = [(float(p[0]), float(p[1])) for p in pts]
        self.S = [0.0]
        for i in range(1, len(self.pts)):
            (x0, y0), (x1, y1) = self.pts[i - 1], self.pts[i]
            self.S.append(self.S[-1] + math.hypot(x1 - x0, y1 - y0))
        xs = [p[0] for p in self.pts]
        ys = [p[1] for p in self.pts]
        self.bbox = (min(xs), min(ys), max(xs), max(ys))

    def project(self, x, y):
        """Return (S, dist, seg) of the closest point on the polyline."""
        best = (None, 1e300, -1)
        for i in range(len(self.pts) - 1):
            ax, ay = self.pts[i]
            bx, by = self.pts[i + 1]
            dx, dy = bx - ax, by - ay
            L2 = dx * dx + dy * dy
            if L2 == 0.0:
                continue
            t = ((x - ax) * dx + (y - ay) * dy) / L2
            t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
            qx, qy = ax + t * dx, ay + t * dy
            d = math.hypot(x - qx, y - qy)
            if d < best[1]:
                best = (self.S[i] + t * math.sqrt(L2), d, i)
        return best

    def tangent(self, seg):
        (ax, ay), (bx, by) = self.pts[seg], self.pts[seg + 1]
        L = math.hypot(bx - ax, by - ay)
        return (bx - ax) / L, (by - ay) / L


def nodes_on_polyline(nodes, poly, tol, S_range=None):
    """nodes: dict label -> (x, y). Returns list of (label, S, seg) sorted by S."""
    x0, y0, x1, y1 = poly.bbox
    out = []
    for lab, (x, y) in nodes.items():
        if x < x0 - tol or x > x1 + tol or y < y0 - tol or y > y1 + tol:
            continue
        S, d, seg = poly.project(x, y)
        if d > tol:
            continue
        if S_range is not None and (S < S_range[0] - tol or S > S_range[1] + tol):
            continue
        out.append((lab, S, seg))
    out.sort(key=lambda r: (r[1], r[0]))
    return out


def node_to_elements(conn, wanted):
    """conn: dict elem -> tuple(node labels). Map wanted node labels to elements."""
    wanted = set(wanted)
    m = {}
    for e, nl in conn.items():
        for n in nl:
            if n in wanted:
                m.setdefault(n, []).append(e)
    return m


def average_element_nodal(values, node_set, elem_set):
    """values: iterable of (elem, node, (sxx, syy, sxy)). Average per node over
    elements belonging to elem_set. Returns dict node -> (sxx, syy, sxy)."""
    acc = {}
    for e, n, s in values:
        if n not in node_set or e not in elem_set:
            continue
        a = acc.setdefault(n, [0.0, 0.0, 0.0, 0])
        a[0] += s[0]; a[1] += s[1]; a[2] += s[2]; a[3] += 1
    return dict((n, (a[0] / a[3], a[1] / a[3], a[2] / a[3])) for n, a in acc.items())


def centroid(conn_e, coords):
    xs = [coords[n][0] for n in conn_e]
    ys = [coords[n][1] for n in conn_e]
    return sum(xs) / len(xs), sum(ys) / len(ys)


def seam_jumps(on_line, coords, u, n2e, side_fn, tol):
    """Pair duplicated seam nodes and return jumps u(+) - u(-).

    on_line : list of (label, S, seg) from nodes_on_polyline
    side_fn : function(label, S, seg) -> +1 (left of the path tangent) or -1
    Returns list of (S, seg, dux, duy) sorted by S. Unpaired nodes (crack tips)
    are skipped; the post-processor adds zero jumps at the tips.
    """
    groups = []
    for lab, S, seg in on_line:
        x, y = coords[lab]
        for g in groups:
            gx, gy = coords[g[0][0]]
            if abs(gx - x) <= tol and abs(gy - y) <= tol:
                g.append((lab, S, seg))
                break
        else:
            groups.append([(lab, S, seg)])
    out = []
    for g in groups:
        if len(g) != 2:
            continue
        (l1, S1, s1), (l2, S2, s2) = g
        sd1, sd2 = side_fn(l1, S1, s1), side_fn(l2, S2, s2)
        if sd1 == sd2:
            continue          # could not resolve sides; skip rather than guess
        lp, lm = (l1, l2) if sd1 > 0 else (l2, l1)
        out.append((0.5 * (S1 + S2), s1, u[lp][0] - u[lm][0], u[lp][1] - u[lm][1]))
    out.sort(key=lambda r: r[0])
    return out


def side_by_centroid(poly, coords, conn, n2e):
    """Side function for straight (bulk) seams: sign of cross(t, c - p)."""
    def f(lab, S, seg):
        tx, ty = poly.tangent(seg)
        px, py = coords[lab]
        s = 0.0
        for e in n2e.get(lab, []):
            cx, cy = centroid(conn[e], coords)
            s += tx * (cy - py) - ty * (cx - px)
        return 1 if s > 0 else -1
    return f


def side_by_material(n2e, elems_plus):
    """Side function for interface seams: + side is material B (left of +S)."""
    def f(lab, S, seg):
        es = n2e.get(lab, [])
        nplus = sum(1 for e in es if e in elems_plus)
        return 1 if nplus * 2 > len(es) else -1
    return f


def sum_nodal_forces(values, node_set, elem_ok):
    """values: iterable of (elem, node, (f1, f2)) from NFORC. Sum per node over the
    elements for which elem_ok(elem) is True (one side of the future crack)."""
    acc = {}
    for e, n, f in values:
        if n not in node_set or not elem_ok(e):
            continue
        a = acc.setdefault(n, [0.0, 0.0])
        a[0] += f[0]; a[1] += f[1]
    return acc


def left_of_line_fn(p0, t, coords, conn):
    """elem_ok for straight lines: element centroid lies left of the tangent t."""
    cache = {}
    def ok(e):
        if e not in cache:
            cx, cy = centroid(conn[e], coords)
            cache[e] = (t[0] * (cy - p0[1]) - t[1] * (cx - p0[0])) > 0.0
        return cache[e]
    return ok
