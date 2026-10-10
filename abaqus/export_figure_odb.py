# -*- coding: utf-8 -*-
"""Export mesh, displacement, stress and (CZM) interface damage of one ODB frame to JSON, for
drawing Abaqus meshes / contour plots in MATLAB. Runs with Abaqus' Python (2017 or 2025):

    abaqus python abaqus/export_figure_odb.py <file.odb> <case.json> <out.json> [t]

t = step time of the frame to export (0..1, nearest available frame; default: last frame).
Output (lengths in the model units):
  nodes: [[label, x, y, ux, uy], ...]
  elements: [[label, material(0 = A, 1 = B), n1, n2, n3, n4, n5, n6], ...]  (CPE6: 3 corner + 3 mid-side)
  stress_nodal: per material, nodal averages of S11, S22, S12 in global axes
  damage: [[x, y, CSDMG], ...] on the cohesive interface (CZM ODBs only)
"""
from __future__ import print_function, division
import sys, json, math
from odbAccess import openOdb
from abaqusConstants import ELEMENT_NODAL

args = [a for a in sys.argv[1:]]
odb_path, case_path, out_path = args[0], args[1], args[2]
t_want = float(args[3]) if len(args) > 3 else None
CASE = json.load(open(case_path))
FE = CASE['fe']


def _d(v):
    try:
        return v.data
    except Exception:
        return v.dataDouble


odb = openOdb(odb_path, readOnly=True)
step = odb.steps[list(odb.steps.keys())[-1]]
frames = step.frames
if t_want is None:
    fr = frames[-1]
else:
    fr = min(frames, key=lambda f: abs(f.frameValue - t_want))
inst = odb.rootAssembly.instances[list(odb.rootAssembly.instances.keys())[0]]
matA = set(e.label for e in inst.elementSets['MAT_A'].elements)
U = dict((v.nodeLabel, _d(v)) for v in fr.fieldOutputs['U'].values)
nodes = [[n.label, float(n.coordinates[0]), float(n.coordinates[1]),
          float(U.get(n.label, (0, 0))[0]), float(U.get(n.label, (0, 0))[1])] for n in inst.nodes]
elements = [[e.label, 0 if e.label in matA else 1] + list(e.connectivity) for e in inst.elements]

out = dict(odb=odb_path, case=CASE.get('name'), frame_time=fr.frameValue, nodes=nodes, elements=elements)

if 'S' in fr.fieldOutputs.keys():
    rot = {}
    for key, mi in (('A', 0), ('B', 1)):
        mk = FE['materials'][key]
        if mk.get('ortho') and abs(mk.get('raster_deg', 0.0)) > 1e-9:
            rot[mi] = math.radians(mk['raster_deg'])
    S = fr.fieldOutputs['S'].getSubset(position=ELEMENT_NODAL)
    lab = list(S.componentLabels)
    i11, i22, i12 = lab.index('S11'), lab.index('S22'), lab.index('S12')
    acc = {}
    for v in S.values:
        mi = 0 if v.elementLabel in matA else 1
        d = _d(v)
        s11, s22, s12 = d[i11], d[i22], d[i12]
        th = rot.get(mi)
        if th is not None:
            c, s_ = math.cos(th), math.sin(th)
            s11, s22, s12 = (c * c * s11 + s_ * s_ * s22 - 2 * c * s_ * s12,
                             s_ * s_ * s11 + c * c * s22 + 2 * c * s_ * s12,
                             c * s_ * (s11 - s22) + (c * c - s_ * s_) * s12)
        a = acc.setdefault((mi, v.nodeLabel), [0.0, 0.0, 0.0, 0])
        a[0] += s11; a[1] += s22; a[2] += s12; a[3] += 1
    out['stress_nodal'] = [[mi, n, a[0] / a[3], a[1] / a[3], a[2] / a[3]] for (mi, n), a in acc.items()]

fk = [k for k in fr.fieldOutputs.keys() if k.split()[0] == 'CSDMG']
if fk:
    coords = dict((n.label, n.coordinates) for n in inst.nodes)
    dmg = []
    for v in fr.fieldOutputs[fk[0]].values:
        x = _d(v)
        x = x[0] if hasattr(x, '__len__') else x
        if v.nodeLabel in coords:
            c = coords[v.nodeLabel]
            dmg.append([float(c[0]), float(c[1]), float(x)])
    out['damage'] = dmg
odb.close()
with open(out_path, 'w') as fh:
    json.dump(out, fh)
print('written', out_path, 'frame time', out['frame_time'])
