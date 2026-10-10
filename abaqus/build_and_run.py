# -*- coding: utf-8 -*-
"""
Abaqus/CAE driver for one FE setup (geometry x material pair x load case).

    abaqus cae noGUI=abaqus/build_and_run.py -- runs/<setup>/case.json [slice_index n_slices]
                                               [uncracked_only] [save_cae]

Workflow
  1. Build the base part: outline, interface polyline and all bulk candidate lines
     are partitioned *with every candidate crack tip as a vertex*. This makes the
     mesh identical in the uncracked and in every cracked model (only the seam
     differs), which keeps Delta(ALLSE) and the crack-closure integral consistent.
  2. Run the uncracked model -> raw/uncracked.json (stresses on all paths, ALLSE, F).
  3. For each path and crack length: copy the base model, assign a seam on the
     crack edges, mesh, run -> raw/<path>__<k>.json (ALLSE, F, crack-face jumps).
  Existing raw files are skipped, so the script can be restarted at any time.

NOTE: written against the Abaqus 2019-2024 scripting API but NOT executed by the
author. Run the uncracked model of one case with "save_cae": true first and check
partitions, seams and the mesh visually before launching a sweep.
"""
from __future__ import print_function, division

import sys
import os
import json
import math
import time
import traceback
import io

from abaqus import mdb, session
from abaqusConstants import *
# In noGUI scripts the Model methods of each CAE module (FieldOutputRequest, DisplacementBC,
# ElemType, Job, ...) are only attached after the module is imported; caeModules imports all.
from caeModules import *
import regionToolset
import mesh
from odbAccess import openOdb

# The star imports above shadow some Python builtins with Abaqus functions (e.g. sum ->
# "arg1; found 'generator'"). Restore the builtins used in this script.
try:
    import __builtin__ as _bi          # Python 2 (Abaqus <= 2023)
except ImportError:
    import builtins as _bi             # Python 3 (Abaqus >= 2024)
for _n in ('sum', 'min', 'max', 'abs', 'round', 'pow', 'any', 'all', 'range', 'len', 'zip', 'map',
           'filter', 'set', 'dict', 'list', 'tuple', 'float', 'int', 'str', 'open', 'sorted'):
    globals()[_n] = getattr(_bi, _n)

# ----------------------------------------------------------------------------
# arguments and imports from the project
# ----------------------------------------------------------------------------
ARGS = [a for a in sys.argv if a.endswith('.json') or a.isdigit()]
CASE_FILE = os.path.abspath([a for a in ARGS if a.endswith('.json')][-1])
NUMS = [int(a) for a in ARGS if a.isdigit()]
SLICE = (NUMS[-2], NUMS[-1]) if len(NUMS) >= 2 else (0, 1)
UNCRACKED_ONLY = 'uncracked_only' in sys.argv      # quick model check (pilot)
# keep=I2f__11,B1A355__05 : re-run only these cracked analyses and KEEP their ODB (figures);
# existing raw results are not touched
KEEP = [t for a in sys.argv if a.startswith('keep=') for t in a[5:].split(',') if t]
FORCE_SAVE_CAE = 'save_cae' in sys.argv            # write uncracked_check.cae for visual inspection

def _to_str(x):
    """Python 2 (Abaqus <= 2023): json gives unicode, several Abaqus API calls require str."""
    if sys.version_info[0] >= 3:
        return x
    if isinstance(x, dict):
        return dict((_to_str(k), _to_str(v)) for k, v in x.items())
    if isinstance(x, list):
        return [_to_str(v) for v in x]
    if isinstance(x, unicode):  # noqa: F821  (py2 only)
        return x.encode('utf-8')
    return x


with open(CASE_FILE) as fh:
    CASE = _to_str(json.load(fh))
# The project may have been moved to another computer after the cases were generated:
# fall back to <project>/<runs dir>/<setup>/case.json -> <project>
PROJECT_ROOT = CASE.get('project_root', '')
if not os.path.isfile(os.path.join(PROJECT_ROOT, 'ffm', 'abq_utils.py')):
    PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(CASE_FILE)))
sys.path.insert(0, PROJECT_ROOT)
from ffm.abq_utils import (Polyline, nodes_on_polyline, node_to_elements, average_element_nodal,
                           seam_jumps, side_by_centroid, side_by_material, sum_nodal_forces,
                           left_of_line_fn)

SETUP_DIR = os.path.dirname(CASE_FILE)
RAW_DIR = os.path.join(SETUP_DIR, 'raw')
WORK_DIR = os.path.join(SETUP_DIR, 'work_s{0}'.format(SLICE[0]))   # one per slice: parallel-safe
for d in (RAW_DIR, WORK_DIR):
    if not os.path.isdir(d):
        os.makedirs(d)
os.chdir(WORK_DIR)

W, H = CASE['W'], CASE['H']
FE = CASE['fe']                    # materials, load, plane, mesh, delta
PERIODIC = CASE.get('lateral', 'free') == 'periodic'
TOL = 1e-6 * max(W, H)
PART, INST = 'SPECIMEN', 'SPECIMEN-1'
IFACE = Polyline(CASE['interface_points'])
BULK_POLYS = [Polyline(bl['points']) for bl in CASE['bulk_lines']]


def log(msg):
    line = '[{0}] {1}'.format(time.strftime('%H:%M:%S'), msg)
    try:
        print(line)
    except Exception:
        pass
    try:
        with open(os.path.join(SETUP_DIR, 'abaqus_log.txt'), 'a') as fh:
            fh.write(line + '\n')
    except UnicodeError:                      # localised OS messages (e.g. Turkish Windows)
        with open(os.path.join(SETUP_DIR, 'abaqus_log.txt'), 'ab') as fh:
            fh.write((line + '\n').encode('ascii', 'replace'))


# ----------------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------------
def _partition_sketch(m, p, polylines):
    try:
        sk = m.ConstrainedSketch(name='__part__', sheetSize=4.0 * max(W, 2 * H))
    except Exception:
        t = p.MakeSketchTransform(sketchPlane=p.faces[0], sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
        sk = m.ConstrainedSketch(name='__part__', sheetSize=4.0 * max(W, 2 * H), transform=t)
    for pl in polylines:
        for (x1, y1), (x2, y2) in zip(pl[:-1], pl[1:]):
            if math.hypot(x2 - x1, y2 - y1) > TOL:
                sk.Line(point1=(x1, y1), point2=(x2, y2))
    p.PartitionFaceBySketch(faces=p.faces, sketch=sk)
    del m.sketches['__part__']


def _ensure_vertex(p, x, y):
    """Make sure a vertex exists at (x, y); partition the edge if the sketch merged it."""
    hit = p.vertices.getClosest(coordinates=((x, y, 0.0),))
    if hit and math.hypot(hit[0][1][0] - x, hit[0][1][1] - y) <= 10 * TOL:
        return
    e = p.edges.getClosest(coordinates=((x, y, 0.0),))[0][0]
    dp = p.DatumPointByCoordinate(coords=(x, y, 0.0))
    p.PartitionEdgeByPoint(edge=e, point=p.datums[dp.id])


def build_base(name='BASE'):
    if name in mdb.models:
        del mdb.models[name]
    m = mdb.Model(name=name)
    s = m.ConstrainedSketch(name='outline', sheetSize=4.0 * max(W, 2 * H))
    s.rectangle(point1=(0.0, -H), point2=(W, H))
    p = m.Part(name=PART, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    p.BaseShell(sketch=s)

    iface_pts = [(b[0], b[1]) for b in CASE['interface_breakpoints']]
    lines = [iface_pts] + [[tuple(q) for q in bl['points']] for bl in CASE['bulk_lines']]
    _partition_sketch(m, p, lines)
    for (x, y) in iface_pts[1:-1]:
        _ensure_vertex(p, x, y)
    for bl in CASE['bulk_lines']:
        for (x, y) in bl['points'][1:]:
            _ensure_vertex(p, x, y)

    # materials and sections (thickness 1 -> energies and forces per unit thickness).
    # Crossing bulk lines (e.g. both neck corners in shear) cut a material region into
    # several faces: classify EVERY face by point-in-polygon against the interface.
    face_groups = _faces_by_material(p)
    for key in ('A', 'B'):
        mat = FE['materials'][key]
        mt = m.Material(name='MAT_' + key)
        oc = mat.get('ortho')
        if oc:
            # axes: 1 = roads, 2 = in-layer transverse, 3 = build direction (out of plane)
            mt.Elastic(type=ENGINEERING_CONSTANTS,
                       table=((oc['E1'], oc['E2'], oc['E3'], oc['nu12'], oc['nu13'], oc['nu23'],
                               oc['G12'], oc['G13'], oc['G23']),))
        else:
            mt.Elastic(table=((mat['E'], mat['nu']),))
        m.HomogeneousSolidSection(name='SEC_' + key, material='MAT_' + key, thickness=1.0)
        p.Set(faces=face_groups[key], name='MAT_' + key)
        p.SectionAssignment(region=p.sets['MAT_' + key], sectionName='SEC_' + key)
        if oc:
            th = math.radians(mat.get('raster_deg', 0.0))
            dc = p.DatumCsysByThreePoints(name='CS_' + key, coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0),
                                          point1=(math.cos(th), math.sin(th), 0.0),
                                          point2=(-math.sin(th), math.cos(th), 0.0))
            p.MaterialOrientation(region=p.sets['MAT_' + key], orientationType=SYSTEM, axis=AXIS_3,
                                  localCsys=p.datums[dc.id], fieldName='',
                                  additionalRotationType=ROTATION_NONE, angle=0.0,
                                  additionalRotationField='', stackDirection=STACK_3)
    p.Set(edges=p.edges.findAt(((0.5 * W, H, 0.0),)), name='TOP')
    p.Set(edges=p.edges.findAt(((0.5 * W, -H, 0.0),)), name='BOT')
    p.Set(vertices=p.vertices.findAt(((0.0, -H, 0.0),)), name='PIN')

    # pre-crack (benchmarks): seam present in every model
    if CASE.get('precrack'):
        _assign_seam(p, CASE['precrack']['seam_points'])

    a = m.rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    inst = a.Instance(name=INST, part=p, dependent=ON)
    m.StaticStep(name='Load', previous='Initial')
    # explicit output requests: default requests are not created in every CAE version/model
    # type, so they are neither relied on nor duplicated
    for rep_name in ('fieldOutputRequests', 'historyOutputRequests'):
        rep = getattr(m, rep_name, None)
        if rep is not None:
            for k in list(rep.keys()):
                del rep[k]
    m.FieldOutputRequest(name='F-FFM', createStepName='Load', variables=('S', 'U', 'RF', 'NFORC'))
    m.HistoryOutputRequest(name='H-ENER', createStepName='Load', variables=('ALLSE', 'ALLWK'))
    d = FE['delta']
    sliding = FE.get('grips', 'clamped') == 'sliding'
    if sliding:
        # frictionless grips (tension only): u2 prescribed, lateral contraction free, one pin
        if FE['load'] != 'tension':
            raise ValueError('sliding grips are only defined for tension')
        m.DisplacementBC(name='BOT', createStepName='Initial', region=inst.sets['BOT'],
                         u1=UNSET, u2=SET, ur3=UNSET)
        m.DisplacementBC(name='PIN', createStepName='Initial', region=inst.sets['PIN'],
                         u1=SET, u2=SET, ur3=UNSET)
        m.DisplacementBC(name='TOP', createStepName='Load', region=inst.sets['TOP'], u1=UNSET, u2=d, ur3=UNSET)
        return m
    m.DisplacementBC(name='BOT', createStepName='Initial', region=inst.sets['BOT'],
                     u1=SET, u2=SET, ur3=UNSET)
    if FE['load'] == 'tension':
        m.DisplacementBC(name='TOP', createStepName='Load', region=inst.sets['TOP'], u1=0.0, u2=d, ur3=UNSET)
    elif FE['load'] == 'shear':
        m.DisplacementBC(name='TOP', createStepName='Load', region=inst.sets['TOP'], u1=d, u2=0.0, ur3=UNSET)
    else:
        raise ValueError(FE['load'])
    return m


def _in_A(x, y):
    """Point-in-polygon: region A = below the interface polyline (closed via the bottom edge)."""
    poly = [tuple(q) for q in CASE['interface_points']] + [(W, -H), (0.0, -H)]
    inside = False
    j = len(poly) - 1
    for i in range(len(poly)):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / (yj - yi) + xi:
            inside = not inside
        j = i
    return inside


def _faces_by_material(p):
    """Face arrays per material built from the face objects themselves (by index): a
    findAt() with the face's own pointOn can fail or hit the neighbour face when that point
    lies on a shared edge (observed for the sine interface)."""
    idx = {'A': [], 'B': []}
    for i, f in enumerate(p.faces):
        x, y = f.pointOn[0][0], f.pointOn[0][1]
        idx['A' if _in_A(x, y) else 'B'].append(i)
    if not idx['A'] or not idx['B']:
        raise RuntimeError('face classification failed: {0} A faces, {1} B faces'
                           .format(len(idx['A']), len(idx['B'])))
    log('faces: {0} in A, {1} in B'.format(len(idx['A']), len(idx['B'])))
    out = {}
    for k, ii in idx.items():
        arr = p.faces[ii[0]:ii[0] + 1]
        for i in ii[1:]:
            arr = arr + p.faces[i:i + 1]
        out[k] = arr
    return out


def _assign_seam(p, seam_points):
    if not seam_points:
        return
    edges = p.edges.findAt(*[((x, y, 0.0),) for x, y in seam_points])
    p.engineeringFeatures.assignSeam(regions=regionToolset.Region(edges=edges))


def _on_lines(x, y):
    if IFACE.project(x, y)[1] <= 10 * TOL:
        return True
    for pl in BULK_POLYS:
        if pl.project(x, y)[1] <= 10 * TOL:
            return True
    return False


def mesh_part(m):
    p = m.parts[PART]
    ms = FE['mesh']
    p.seedPart(size=ms['h_global'], deviationFactor=0.1, minSizeFactor=0.1)
    # every short edge on interface / bulk lines gets n_per_edge elements -> geometric grading
    for e in p.edges:
        x, y, _ = e.pointOn[0]
        if not _on_lines(x, y):
            continue
        L = e.getSize(printResults=False)
        n = max(ms['n_per_edge'], int(math.ceil(L / ms['h_line'])))
        p.seedEdgeByNumber(edges=(e,), number=n, constraint=FINER)
    if PERIODIC:
        # identical seeds on the left and right cell edges -> matching node pairs
        for e in p.edges:
            x, y, _ = e.pointOn[0]
            if abs(x) < TOL or abs(x - W) < TOL:
                L = e.getSize(printResults=False)
                # round before ceil: lengths that are exact multiples of h_global (e.g. 8.5/0.5)
                # come out as 16.9999999 on one side and 17.0000001 on the other -> one element
                # more on one cell edge and no matching node pairs (seen for sine_ar1)
                n = int(math.ceil(round(L / ms['h_global'], 6)))
                p.seedEdgeByNumber(edges=(e,), number=max(4, n), constraint=FIXED)
    p.setMeshControls(regions=p.faces, elemShape=TRI, technique=FREE)
    if FE['plane'] == 'strain':
        et = (mesh.ElemType(elemCode=CPE8, elemLibrary=STANDARD), mesh.ElemType(elemCode=CPE6, elemLibrary=STANDARD))
    else:
        et = (mesh.ElemType(elemCode=CPS8, elemLibrary=STANDARD), mesh.ElemType(elemCode=CPS6, elemLibrary=STANDARD))
    p.setElementType(regions=(p.faces,), elemTypes=et)
    p.generateMesh()
    m.rootAssembly.regenerate()
    if PERIODIC:
        add_periodic(m)
    return len(p.elements)


def add_periodic(m):
    """Interior cell of a periodic array: u(W, y) - u(0, y) = u_RP for every node pair on the
    left/right cell edges (top/bottom corners excluded: they carry the grip BCs). The dummy
    reference point RP has free dofs = macroscopic stretch / shear of the array across the
    cell, so no lateral force is imposed (sigma_xx, sigma_xy macro = 0 on the cell sides).
    A crack that reaches a cell edge duplicates the boundary node there (seam): every
    duplicate is tied to the single partner node, i.e. the crack ends at the cell edge."""
    p = m.parts[PART]
    a = m.rootAssembly
    tol = 1e-6 * max(W, H)
    key = lambda y: int(round(y / (1e-5 * H)))

    def side(x0):
        d = {}
        for n in p.nodes:
            x, y = n.coordinates[0], n.coordinates[1]
            if abs(x - x0) < tol and abs(abs(y) - H) > tol:
                d.setdefault(key(y), []).append(n.label)
        return d
    left, right = side(0.0), side(W)
    if sorted(left) != sorted(right):
        raise RuntimeError('periodic BC: left/right cell-edge nodes do not match ({0} vs {1} positions)'
                           .format(len(left), len(right)))
    rp = a.ReferencePoint(point=(1.2 * W, 0.0, 0.0))
    a.Set(referencePoints=(a.referencePoints[rp.id],), name='PBC_RP')
    eqs = []                                    # (dependent set, independent set, sign of RP)
    # CZM: the cohesive seam crosses BOTH cell edges, i.e. it continues into the neighbouring
    # cells. Its duplicated nodes must then be paired face by face (A with A, B with B);
    # tying them all together would pin the two faces at the cell edge.
    by_face = bool(CASE.get('czm'))
    if by_face:
        nodesA = set()
        for e in p.sets['MAT_A'].elements:
            nodesA.update(nd.label for nd in e.getNodes())
    for k in sorted(left):
        L, R = left[k], right[k]
        if by_face and len(L) == len(R) > 1:
            L = sorted(L, key=lambda n: n not in nodesA)          # A-face node first
            R = sorted(R, key=lambda n: n not in nodesA)
            for nl, nr in zip(L, R):
                eqs.append((('R', nr), ('L', nl), -1.0))          # u_R = u_L + u_RP, per face
            continue
        eqs.append((('R', R[0]), ('L', L[0]), -1.0))       # u_R0 = u_L0 + u_RP
        for nl in L[1:]:
            eqs.append((('L', nl), ('L', L[0]), 0.0))       # u_L = u_L0 (eliminated dofs never reused)
        for nr in R[1:]:
            eqs.append((('R', nr), ('L', L[0]), -1.0))      # u_R = u_L0 + u_RP
    made = set()

    def nset(tag, label):
        nm = 'P{0}{1}'.format(tag, label)
        if nm not in made:
            a.SetFromNodeLabels(name=nm, nodeLabels=((INST, (label,)),))
            made.add(nm)
        return nm
    for i, (dep, ind, sgn) in enumerate(eqs):
        sd, si = nset(*dep), nset(*ind)
        for dof in (1, 2):
            terms = ((1.0, sd, dof), (-1.0, si, dof))
            if sgn:
                terms += ((sgn, 'PBC_RP', dof),)
            m.Equation(name='PBC{0}_{1}'.format(i, dof), terms=terms)
    log('periodic cell: {0} node pairs tied ({1} equations per dof)'.format(len(left), len(eqs)))


def _crack_geometry(path, crack):
    """Function (x, y) -> True if the point lies on THIS crack (used to reject duplicated
    nodes that do not belong to the crack faces)."""
    tol = 1e-4 * max(W, H) / 8.0
    if path['kind'] == 'bulk':
        ox, oy = path['origin']
        ex, ey = math.cos(math.radians(path['angle_deg'])), math.sin(math.radians(path['angle_deg']))
        da = crack['da']

        def on(x, y):
            s_ = (x - ox) * ex + (y - oy) * ey
            d = abs((y - oy) * ex - (x - ox) * ey)
            return d < tol and -tol <= s_ <= da + tol
        return on
    poly = Polyline(CASE['interface_points'])
    lo, hi = crack['S_interval']

    def on(x, y):
        S, dist, seg = poly.project(x, y)
        return dist < tol and lo - tol <= S <= hi + tol
    return on


def add_crack_contact(m, path, crack, mode='hard'):
    """Frictionless contact between the two crack faces (zero initial gap), so closing faces
    cannot interpenetrate. Frictionless contact with zero gap keeps the solution positively
    homogeneous in the load (lambda scaling of the post-processing stays valid).
    Crack-face element edges: quadratic seam edges have a duplicated mid-side node lying on
    this crack. Side: material (interface crack) or side of the crack line (bulk crack), judged
    locally from the element centroid relative to its own crack-face edge.
    mode='hard' (direct, default) or 'penalty' (retry when the hard contact does not converge).
    Returns False (contact skipped) if no valid two-sided surface can be built."""
    p = m.parts[PART]
    tol = 1e-7 * max(W, H)
    pos = {}
    for n in p.nodes:
        k = (int(round(n.coordinates[0] / tol)), int(round(n.coordinates[1] / tol)))
        pos.setdefault(k, []).append(n)
    on_crack = _crack_geometry(path, crack)
    dup = set(nd.label for v in pos.values() if len(v) > 1 for nd in v
              if on_crack(nd.coordinates[0], nd.coordinates[1]))
    if not dup:
        log('crack contact: no seam nodes found - skipped')
        return False
    matA = set(e.label for e in p.sets['MAT_A'].elements)
    if path['kind'] == 'bulk':
        ex, ey = math.cos(math.radians(path['angle_deg'])), math.sin(math.radians(path['angle_deg']))
    faces = {1: {1: [], 2: [], 3: []}, 2: {1: [], 2: [], 3: []}}
    use = {1: {}, 2: {}}                       # node -> number of faces using it, per side
    for e in p.elements:
        nodes = e.getNodes()
        if len(nodes) < 6:
            continue
        for fi, (ia, ib, mid) in ((1, (0, 1, nodes[3])), (2, (1, 2, nodes[4])), (3, (2, 0, nodes[5]))):
            if mid.label not in dup:
                continue
            if path['kind'] == 'interface':
                side = 1 if e.label in matA else 2
            else:
                cx = sum(nd.coordinates[0] for nd in nodes[:3]) / 3.0
                cy = sum(nd.coordinates[1] for nd in nodes[:3]) / 3.0
                mx, my = mid.coordinates[0], mid.coordinates[1]
                side = 1 if ex * (cy - my) - ey * (cx - mx) > 0 else 2
            faces[side][fi].append(e.label)
            for nd in (nodes[ia], nodes[ib], mid):
                use[side][nd.label] = use[side].get(nd.label, 0) + 1
    n1 = sum(len(v) for v in faces[1].values())
    n2 = sum(len(v) for v in faces[2].values())
    bad = [(sd, l) for sd in (1, 2) for l, c in use[sd].items() if c > 2]
    if not n1 or not n2 or bad:
        log('crack contact: invalid crack-face surface ({0}/{1} edges, {2} non-manifold nodes) - '
            'contact skipped for this analysis'.format(n1, n2, len(bad)))
        return False
    for side in (1, 2):
        if 'CF{0}'.format(side) in p.surfaces.keys():
            continue                                  # retry: surfaces already built
        kw = {}
        for fi in (1, 2, 3):
            if faces[side][fi]:
                kw['face{0}Elements'.format(fi)] = p.elements.sequenceFromLabels(faces[side][fi])
        p.Surface(name='CF{0}'.format(side), **kw)
    a = m.rootAssembly
    a.regenerate()
    inst = a.instances[INST]
    pname = 'CF_PROP_' + mode.upper()
    if mode == 'cohesive' and pname not in m.interactionProperties.keys():
        cz = CASE['czm']
        pr = m.ContactProperty(pname)
        K = cz['K']
        pr.CohesiveBehavior(defaultPenalties=OFF, table=((K, K, K),))
        law = cz.get('mixed_mode', 'power')
        dmg = dict(criterion=QUAD_TRACTION, initTable=((cz['sigma_c'], cz['tau_c'], cz['tau_c']),),
                   useEvolution=ON, evolutionType=ENERGY, softening=LINEAR, useMixedMode=ON,
                   modeMixRatio=ENERGY, evolTable=((cz['G_Ic'], cz['G_IIc'], cz['G_IIc']),),
                   useStabilization=ON, viscosityCoef=cz.get('viscosity', 1e-5))
        if law == 'bk':
            dmg.update(mixedModeType=BK, exponent=cz.get('bk_eta', 2.0))
        else:                                   # power law, exponent 1 = linear (as in the FFM)
            dmg.update(mixedModeType=POWER_LAW, exponent=cz.get('power', 1.0))
        pr.Damage(**dmg)
        pr.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON)   # after full damage
        pr.TangentialBehavior(formulation=FRICTIONLESS)
    if pname not in m.interactionProperties.keys():
        pr = m.ContactProperty(pname)
        if mode == 'penalty':
            pr.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON,
                              constraintEnforcementMethod=PENALTY)
        else:
            pr.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON)
        pr.TangentialBehavior(formulation=FRICTIONLESS)
    common = dict(name='CRACK_FACES', createStepName='Initial', sliding=SMALL,
                  interactionProperty=pname, adjustMethod=NONE, initialClearance=OMIT,
                  thickness=OFF)
    try:        # Abaqus 2022+ : main / secondary
        m.SurfaceToSurfaceContactStd(main=inst.surfaces['CF1'], secondary=inst.surfaces['CF2'], **common)
    except TypeError:   # Abaqus <= 2021 : master / slave
        m.SurfaceToSurfaceContactStd(master=inst.surfaces['CF1'], slave=inst.surfaces['CF2'], **common)
    if mode == 'cohesive':
        cz = CASE['czm']
        m.steps['Load'].setValues(initialInc=cz.get('inc0', 0.005), minInc=cz.get('inc_min', 1e-9),
                                  maxInc=cz.get('inc_max', 0.005), maxNumInc=cz.get('n_inc_max', 20000))
        try:
            m.steps['Load'].control.setValues(allowPropagation=OFF, resetDefaultValues=OFF,
                                              discontinuous=ON)
        except Exception:
            pass
    elif mode == 'penalty':
        m.steps['Load'].setValues(initialInc=0.05, minInc=1e-10, maxInc=0.25, maxNumInc=2000)
        try:
            m.steps['Load'].control.setValues(allowPropagation=OFF, resetDefaultValues=OFF,
                                              discontinuous=ON)
        except Exception:
            pass
    else:
        m.steps['Load'].setValues(initialInc=0.5, minInc=1e-6, maxInc=1.0, maxNumInc=100)
    log('crack contact ({0}): {1} + {2} crack-face edges'.format(mode, n1, n2))
    return True


def run_job(m, jobname):
    for ext in ('.lck',):
        if os.path.exists(jobname + ext):
            os.remove(jobname + ext)
    j = mdb.Job(name=jobname, model=m.name, numCpus=FE.get('ncpus', 1), memory=90, memoryUnits=PERCENTAGE,
                nodalOutputPrecision=FULL)   # double-precision U/RF: small crack openings
    j.submit(consistencyChecking=OFF)
    j.waitForCompletion()
    del mdb.jobs[jobname]
    _check_solver(jobname)
    return os.path.join(WORK_DIR, jobname + '.odb')


def _check_solver(jobname):
    """Raise with the solver's own ***ERROR lines if Abaqus/Standard did not finish."""
    msgs = []
    for ext in ('.dat', '.msg'):
        fn = os.path.join(WORK_DIR, jobname + ext)
        if os.path.exists(fn):
            with io.open(fn, 'r', errors='ignore') as fh:     # py2 + py3, any code page
                lines = fh.readlines()
            for i, ln in enumerate(lines):
                if '***ERROR' in ln:
                    msgs.append(ext + ': ' + ' '.join(x.strip() for x in lines[i:i + 4]))
    sta = os.path.join(WORK_DIR, jobname + '.sta')
    done = False
    if os.path.exists(sta):
        with io.open(sta, 'r', errors='ignore') as fh:
            done = 'COMPLETED SUCCESSFULLY' in fh.read()
    if msgs or not done:
        raise RuntimeError('Abaqus/Standard failed for {0}:\n  '.format(jobname) +
                           ('\n  '.join(msgs[:8]) if msgs else 'no ***ERROR found; see .dat/.msg/.sta'))


# ----------------------------------------------------------------------------
# extraction
# ----------------------------------------------------------------------------
def _read_basic(odb):
    step = odb.steps['Load']
    frame = step.frames[-1]
    inst = odb.rootAssembly.instances[INST]
    coords = dict((n.label, (n.coordinates[0], n.coordinates[1])) for n in inst.nodes)
    conn = dict((e.label, tuple(e.connectivity)) for e in inst.elements)
    elemA = set(e.label for e in inst.elementSets['MAT_A'].elements)
    elemB = set(e.label for e in inst.elementSets['MAT_B'].elements)
    hr = [r for k, r in step.historyRegions.items() if k.startswith('Assembly')][0]
    allse = hr.historyOutputs['ALLSE'].data[-1][1]
    allwk = hr.historyOutputs['ALLWK'].data[-1][1]
    rf = frame.fieldOutputs['RF'].getSubset(region=inst.nodeSets['TOP'])
    comp = 1 if FE['load'] == 'tension' else 0
    F = sum(_d(v)[comp] for v in rf.values)
    return frame, inst, coords, conn, elemA, elemB, allse, allwk, F


def _d(v):
    """FieldValue data in either precision. With nodalOutputPrecision=FULL the nodal fields
    (U, RF, possibly NFORC) are double precision and v.data raises 'Use dataDouble()'."""
    try:
        return v.data
    except Exception:
        return v.dataDouble


def _plain(x):
    """ODB values are numpy/Abaqus scalar types (float32 ...) that json cannot write."""
    if isinstance(x, dict):
        return dict((str(k), _plain(v)) for k, v in x.items())
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if x is None or isinstance(x, (bool, str)):
        return x
    if isinstance(x, int) or type(x).__name__ in ('long', 'int32', 'int64', 'uint32'):
        return int(x)
    try:
        return float(x)
    except (TypeError, ValueError):
        return x


def _is_done(path):
    """A result counts as finished only if it is valid JSON (guards against truncated files)."""
    if not os.path.exists(path):
        return False
    try:
        with open(path) as fh:
            json.load(fh)
        return True
    except ValueError:
        log('invalid result file, recomputing: ' + os.path.basename(path))
        os.remove(path)
        return False


def _write_json(res, out_path):
    """Atomic write: a crash while writing never leaves a truncated result that a restart
    would mistake for a finished analysis."""
    tmp = out_path + '.tmp'
    with open(tmp, 'w') as fh:
        json.dump(_plain(res), fh)
    if os.path.exists(out_path):
        os.remove(out_path)
    os.rename(tmp, out_path)


def _read_nforc(frame, want):
    """Element-nodal forces [(elem, node, (F1, F2))] for nodes in `want`.
    The ODB stores NFORC either as one vector field 'NFORC' or as scalar fields
    'NFORC1', 'NFORC2' (common); both layouts are handled."""
    fo = frame.fieldOutputs
    keys = list(fo.keys())
    if 'NFORC1' in keys and 'NFORC2' in keys:
        f1 = dict(((v.elementLabel, v.nodeLabel), _d(v)) for v in fo['NFORC1'].values if v.nodeLabel in want)
        f2 = dict(((v.elementLabel, v.nodeLabel), _d(v)) for v in fo['NFORC2'].values if v.nodeLabel in want)
        return [(k[0], k[1], (f1[k], f2.get(k, 0.0))) for k in f1]
    if 'NFORC' in keys:
        NF = fo['NFORC']
        nlab = list(NF.componentLabels)
        j1, j2 = nlab.index('NFORC1'), nlab.index('NFORC2')
        return [(v.elementLabel, v.nodeLabel, (_d(v)[j1], _d(v)[j2])) for v in NF.values if v.nodeLabel in want]
    raise KeyError('NFORC not found in ODB field outputs: {0}'.format(keys))


def extract_uncracked(odb_path, out_path):
    odb = openOdb(odb_path, readOnly=True)
    frame, inst, coords, conn, elemA, elemB, allse, allwk, F = _read_basic(odb)
    S = frame.fieldOutputs['S'].getSubset(position=ELEMENT_NODAL)
    labels = list(S.componentLabels)          # plane strain: S11,S22,S33,S12; plane stress: S11,S22,S12
    i11, i22, i12 = labels.index('S11'), labels.index('S22'), labels.index('S12')
    want = {}
    on_if = nodes_on_polyline(coords, IFACE, 10 * TOL)
    for lab, s_, seg in on_if:
        want[lab] = 1
    bulk = {}
    for bl in CASE['bulk_lines']:
        pl = Polyline(bl['points'])
        bulk[bl['id']] = (pl, nodes_on_polyline(coords, pl, 10 * TOL), bl['material'])
        for lab, s_, seg in bulk[bl['id']][1]:
            want[lab] = 1
    # with a material orientation Abaqus reports S in the material axes -> rotate to global
    rot = {}
    for key, eset in (('A', elemA), ('B', elemB)):
        mk = FE['materials'][key]
        if mk.get('ortho') and abs(mk.get('raster_deg', 0.0)) > 1e-9:
            rot[key] = math.radians(mk['raster_deg'])

    def to_global(e, s11, s22, s12):
        th = rot.get('A' if e in elemA else 'B')
        if th is None:
            return (s11, s22, s12)
        c, s_ = math.cos(th), math.sin(th)
        return (c * c * s11 + s_ * s_ * s22 - 2 * c * s_ * s12,
                s_ * s_ * s11 + c * c * s22 + 2 * c * s_ * s12,
                c * s_ * (s11 - s22) + (c * c - s_ * s_) * s12)
    vals = [(v.elementLabel, v.nodeLabel, to_global(v.elementLabel, *[_d(v)[i] for i in (i11, i22, i12)]))
            for v in S.values if v.nodeLabel in want]
    node_set = set(want.keys())
    sA = average_element_nodal(vals, node_set, elemA)
    sB = average_element_nodal(vals, node_set, elemB)
    # nodal forces transmitted across the future crack line (for discrete crack closure):
    # sum of NFORC over the elements on the + side (B for the interface, left side for bulk lines)
    fvals = _read_nforc(frame, want)
    fB = sum_nodal_forces(fvals, node_set, lambda e: e in elemB)
    res = dict(ALLSE=allse, ALLWK=allwk, F=F, n_elem=len(conn), interface=[], bulk={})
    for lab, s_, seg in on_if:
        res['interface'].append([s_, seg, coords[lab][0], coords[lab][1],
                                 list(sA.get(lab, (None,) * 3)), list(sB.get(lab, (None,) * 3)),
                                 list(fB.get(lab, (0.0, 0.0)))])
    for bid, (pl, on, mat) in bulk.items():
        src = sA if mat == 'A' else sB
        elems_m = elemA if mat == 'A' else elemB
        t = pl.tangent(0)
        left = left_of_line_fn(pl.pts[0], t, coords, conn)
        fL = sum_nodal_forces(fvals, set(r[0] for r in on), lambda e: (e in elems_m) and left(e))
        res['bulk'][bid] = [[s_, seg] + list(src[lab]) + list(fL.get(lab, (0.0, 0.0)))
                            for lab, s_, seg in on if lab in src]
    odb.close()
    _write_json(res, out_path)


def extract_cracked(odb_path, out_path, path, crack):
    odb = openOdb(odb_path, readOnly=True)
    frame, inst, coords, conn, elemA, elemB, allse, allwk, F = _read_basic(odb)
    u = dict((v.nodeLabel, tuple(_d(v)[:2])) for v in frame.fieldOutputs['U'].values)
    if path['kind'] == 'interface':
        poly = IFACE
        on = nodes_on_polyline(coords, poly, 10 * TOL, S_range=crack['S_interval'])
        n2e = node_to_elements(conn, [r[0] for r in on])
        side = side_by_material(n2e, elemB)
    else:
        poly = Polyline(crack['line'])
        on = nodes_on_polyline(coords, poly, 10 * TOL)
        n2e = node_to_elements(conn, [r[0] for r in on])
        side = side_by_centroid(poly, coords, conn, n2e)
    jumps = seam_jumps(on, coords, u, n2e, side, 10 * TOL)
    res = dict(ALLSE=allse, ALLWK=allwk, F=F, n_elem=len(conn), jumps=[list(j) for j in jumps])
    odb.close()
    _write_json(res, out_path)


def _cleanup(jobname, keep_odb):
    for ext in ('.com', '.dat', '.msg', '.prt', '.sim', '.sta', '.log', '.ipm', '.inp', '.lck'):
        f = jobname + ext
        if os.path.exists(f):
            try:
                os.remove(f)
            except OSError:
                pass
    if not keep_odb and os.path.exists(jobname + '.odb'):
        os.remove(jobname + '.odb')


# ----------------------------------------------------------------------------
# main loop
# ----------------------------------------------------------------------------
def main_logged():
    try:
        main()
    except Exception:
        log('FATAL ' + traceback.format_exc().replace('\n', ' | '))
        raise


def czm_main():
    """Cohesive-zone counterpart of a configuration: the whole interface is a seam whose two
    faces are joined by a surface-based cohesive law (quadratic traction initiation with the
    interface strengths, linear softening with the interface toughnesses). Displacement-
    controlled to CASE['fe']['delta']; the load-displacement curve and the interface damage are
    written to raw/czm.json (also when the solver stops early, e.g. at the final collapse)."""
    out = os.path.join(RAW_DIR, 'czm.json')
    if _is_done(out):
        log('czm: done')
        return
    if SLICE[0] != 0:                        # one CZM analysis per case: other slices idle
        log('czm: slice {0} has nothing to do'.format(SLICE[0]))
        return
    log('czm case {0}'.format(CASE['name']))
    m = build_base('CZM')
    # load-displacement at every increment (history, small); damage fields every 10th increment
    m.fieldOutputRequests['F-FFM'].setValues(variables=('U', 'RF', 'CSDMG', 'CSQUADSCRT'), frequency=10)
    m.historyOutputRequests['H-ENER'].setValues(variables=('ALLSE', 'ALLWK', 'ALLVD', 'ALLSD', 'ALLDMD', 'ALLCD'))
    m.HistoryOutputRequest(name='H-TOP', createStepName='Load', variables=('U1', 'U2', 'RF1', 'RF2'),
                           region=m.rootAssembly.instances[INST].sets['TOP'])
    bp = [tuple(b[:2]) for b in CASE['interface_breakpoints']]
    seam = [(0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])) for a, b in zip(bp[:-1], bp[1:])]
    _assign_seam(m.parts[PART], seam)
    ne = mesh_part(m)
    ip = CASE['interface_points']
    S_tot = sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(ip[:-1], ip[1:]))
    whole = dict(kind='interface', mode='edge', direction=1)
    ok = add_crack_contact(m, whole, dict(S_interval=[-1e-9, S_tot + 1e-9], da=S_tot), mode='cohesive')
    if not ok:
        raise RuntimeError('cohesive surfaces could not be built')
    log('czm model: {0} elements'.format(ne))
    if FE.get('save_cae') or FORCE_SAVE_CAE:
        mdb.saveAs(pathName=os.path.join(SETUP_DIR, 'czm_check.cae'))
    completed = True
    try:
        odb_path = run_job(m, 'czm')
    except RuntimeError as err:
        completed = False                   # keep what was computed (peak usually reached)
        odb_path = os.path.join(WORK_DIR, 'czm.odb')
        log('czm: solver stopped early - extracting the computed part ({0})'.format(str(err).splitlines()[-1][:120]))
        if not os.path.exists(odb_path):
            raise
    extract_czm(odb_path, out, completed)
    _cleanup('czm', keep_odb=True)


def extract_czm(odb_path, out_path, completed):
    odb = openOdb(odb_path, readOnly=True)
    step = odb.steps['Load']
    comp = 'RF2' if FE['load'] == 'tension' else 'RF1'
    ucomp = 'U2' if FE['load'] == 'tension' else 'U1'
    # load-displacement from the node history of the loaded edge (every increment)
    Fsum, U = {}, {}
    for k, hr in step.historyRegions.items():
        if not k.startswith('Node'):
            continue
        ho = hr.historyOutputs
        if comp in ho.keys():
            for t, v in ho[comp].data:
                Fsum[t] = Fsum.get(t, 0.0) + v
        if ucomp in ho.keys() and not U:
            U = dict((t, v) for t, v in ho[ucomp].data)
    times = sorted(Fsum)
    # interface damage from the field frames (every 10th increment)
    dmg = []
    for fr in step.frames:
        row = dict(t=fr.frameValue)
        for key, name in (('CSDMG', 'dmg'), ('CSQUADSCRT', 'init')):
            # contact outputs are stored per contact pair: 'CSDMG    ASSEMBLY_..._CF2/...'
            fk = [k for k in fr.fieldOutputs.keys() if k.split()[0] == key]
            if fk:
                vals = []
                for v in fr.fieldOutputs[fk[0]].values:
                    x = _d(v)
                    vals.append(x[0] if hasattr(x, '__len__') else x)
                if vals:
                    row[name + '_max'] = max(vals)
                    if key == 'CSDMG':
                        row['frac_damaged'] = sum(1 for x in vals if x > 1e-6) / float(len(vals))
                        row['frac_failed'] = sum(1 for x in vals if x >= 0.99) / float(len(vals))
        dmg.append(row)
    frames = []
    j = 0
    for t in times:
        while j + 1 < len(dmg) and dmg[j + 1]['t'] <= t + 1e-12:
            j += 1
        row = dict(t=t, u=U.get(t, 0.0), F=Fsum[t])
        if dmg:
            for k2 in ('frac_failed', 'frac_damaged', 'dmg_max', 'init_max'):
                if k2 in dmg[j] and dmg[j]['t'] <= t + 1e-12:
                    row[k2] = dmg[j][k2]
        frames.append(row)
    hr = [r for k, r in step.historyRegions.items() if k.startswith('Assembly')][0]
    hist = {}
    for k in ('ALLSE', 'ALLWK', 'ALLVD', 'ALLSD', 'ALLDMD', 'ALLCD'):
        if k in hr.historyOutputs.keys():
            hist[k] = [list(x) for x in hr.historyOutputs[k].data]
    odb.close()
    _write_json(dict(name=CASE['name'], completed=completed, W=W, load=FE['load'],
                     frames=frames, damage_frames=dmg, history=hist, czm=CASE['czm']), out_path)
    log('czm: {0} increments, peak F {1:.4g}, completed={2}'.format(
        len(frames), max(abs(f['F']) for f in frames) if frames else 0.0, completed))


def main():
    if CASE.get('czm'):
        return czm_main()
    log('case {0}  slice {1}/{2}'.format(CASE['name'], SLICE[0], SLICE[1]))
    base = build_base()
    # uncracked (only slice 0 runs it)
    out0 = os.path.join(RAW_DIR, 'uncracked.json')
    if SLICE[0] == 0 and (not _is_done(out0) or 'uncracked' in KEEP):
        m = mdb.Model(name='UNCR', objectToCopy=base)
        ne = mesh_part(m)
        log('uncracked: {0} elements'.format(ne))
        if FE.get('save_cae') or FORCE_SAVE_CAE:
            # a .cae that is open in Abaqus/CAE is locked: fall back to a new name, never abort
            for nm in ('uncracked_check.cae', 'uncracked_check_{0}.cae'.format(time.strftime('%H%M%S'))):
                try:
                    mdb.saveAs(pathName=os.path.join(SETUP_DIR, nm))
                    log('model saved: ' + nm)
                    break
                except Exception as e:
                    log('could not save {0} ({1}); trying another name'.format(nm, str(e).splitlines()[0]))
        odb = run_job(m, 'uncr')
        if not _is_done(out0):
            extract_uncracked(odb, out0)
        _cleanup('uncr', keep_odb=True)
        del mdb.models['UNCR']

    if UNCRACKED_ONLY:
        log('uncracked_only: done')
        return
    jobs = []
    for path in CASE['paths']:
        for k, crack in enumerate(path['cracks']):
            jobs.append((path, k, crack))
    jobs = [j for i, j in enumerate(jobs) if i % SLICE[1] == SLICE[0]]
    for path, k, crack in jobs:
        tag = '{0}__{1:02d}'.format(path['pid'], k)
        out = os.path.join(RAW_DIR, tag + '.json')
        if KEEP:
            if tag not in KEEP:
                continue
        elif _is_done(out):
            continue
        t0 = time.time()
        try:
            mname = 'C_' + tag
            m = mdb.Model(name=mname, objectToCopy=base)
            _assign_seam(m.parts[PART], list(crack['seam_points']))   # pre-crack seam is already in BASE
            mesh_part(m)
            contact = bool(FE.get('crack_contact')) and add_crack_contact(m, path, crack)
            jname = 'c_' + tag.replace('.', '_')
            try:
                odb = run_job(m, jname)
            except RuntimeError as err:
                if not (contact and 'TOO MANY ATTEMPTS' in str(err)):
                    raise
                # hard contact did not converge: same crack again with penalty contact,
                # small increments and the discontinuous-analysis controls
                log('{0}: hard contact did not converge - retry with penalty contact'.format(tag))
                del m.interactions['CRACK_FACES']
                add_crack_contact(m, path, crack, mode='penalty')
                jname = jname + '_pen'
                odb = run_job(m, jname)
            if not (KEEP and _is_done(out)):
                extract_cracked(odb, out, path, crack)
            _cleanup(jname, keep_odb=bool(KEEP) or FE.get('keep_odb', False))
            if KEEP:
                log('{0}: ODB kept -> {1}'.format(tag, os.path.join(WORK_DIR, jname + '.odb')))
            del mdb.models[mname]
            log('{0} da={1:.4g} ok ({2:.1f}s)'.format(tag, crack['da'], time.time() - t0))
        except Exception:
            log('{0} FAILED\n{1}'.format(tag, traceback.format_exc()))
            try:
                del mdb.models['C_' + tag]
            except Exception:
                pass


main_logged()
