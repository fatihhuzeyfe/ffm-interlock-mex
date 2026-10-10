# -*- coding: utf-8 -*-
"""Export the stress/displacement field of one ODB to a compact JSON file for plotting with
matplotlib (scripts/plot_fields.py). Runs with Abaqus' Python (odbAccess), 2017 or 2025:

    abaqus python abaqus/export_field.py <file.odb> <case.json> <out.json>

Stresses are returned in GLOBAL axes (material-orientation output is rotated back with the
road angle of each material) and averaged at the nodes separately in each material, so the
stress jump across the bimaterial interface is preserved. 6-node triangles are exported with
their corner and mid-side nodes; the plotting script splits them into 4 linear triangles.
"""
from __future__ import print_function, division
import sys, json, math
from odbAccess import openOdb
from abaqusConstants import ELEMENT_NODAL

odb_path, case_path, out_path = sys.argv[-3:]
CASE = json.load(open(case_path))
FE = CASE['fe']


def _d(v):
    try:
        return v.data
    except Exception:
        return v.dataDouble


odb = openOdb(odb_path, readOnly=True)
step = odb.steps[list(odb.steps.keys())[-1]]
frame = step.frames[-1]
inst = odb.rootAssembly.instances[list(odb.rootAssembly.instances.keys())[0]]
coords = dict((n.label, (float(n.coordinates[0]), float(n.coordinates[1]))) for n in inst.nodes)
conn = dict((e.label, tuple(e.connectivity)) for e in inst.elements)
mat_of = {}
for key in ('A', 'B'):
    for e in inst.elementSets['MAT_' + key].elements:
        mat_of[e.label] = key

rot = {}
for key in ('A', 'B'):
    mk = FE['materials'][key]
    if mk.get('ortho') and abs(mk.get('raster_deg', 0.0)) > 1e-9:
        rot[key] = math.radians(mk['raster_deg'])

S = frame.fieldOutputs['S'].getSubset(position=ELEMENT_NODAL)
lab = list(S.componentLabels)
i11, i22, i12 = lab.index('S11'), lab.index('S22'), lab.index('S12')
i33 = lab.index('S33') if 'S33' in lab else None
acc = {}
for v in S.values:
    m = mat_of.get(v.elementLabel, 'A')
    d = _d(v)
    s11, s22, s12 = d[i11], d[i22], d[i12]
    th = rot.get(m)
    if th is not None:
        c, s_ = math.cos(th), math.sin(th)
        s11, s22, s12 = (c * c * s11 + s_ * s_ * s22 - 2 * c * s_ * s12,
                         s_ * s_ * s11 + c * c * s22 + 2 * c * s_ * s12,
                         c * s_ * (s11 - s22) + (c * c - s_ * s_) * s12)
    s33 = d[i33] if i33 is not None else 0.0
    a = acc.setdefault((m, v.nodeLabel), [0.0, 0.0, 0.0, 0.0, 0])
    a[0] += s11; a[1] += s22; a[2] += s12; a[3] += s33; a[4] += 1
U = dict((v.nodeLabel, tuple(float(x) for x in _d(v)[:2])) for v in frame.fieldOutputs['U'].values)

out = {'case': CASE.get('name'), 'odb': odb_path, 'materials': {}}
for key in ('A', 'B'):
    els = [e for e in conn if mat_of.get(e) == key]
    nodes = sorted(set(n for e in els for n in conn[e]))
    idx = dict((n, i) for i, n in enumerate(nodes))
    stress = []
    for n in nodes:
        a = acc.get((key, n), [0, 0, 0, 0, 1])
        stress.append([a[0] / a[4], a[1] / a[4], a[2] / a[4], a[3] / a[4]])
    out['materials'][key] = {
        'name': FE['materials'][key]['name'],
        'xy': [coords[n] for n in nodes],
        'u': [U.get(n, (0.0, 0.0)) for n in nodes],
        'stress': stress,                                   # S11, S22, S12, S33 (global axes)
        'tri6': [[idx[n] for n in conn[e]] for e in els],
    }
odb.close()
with open(out_path, 'w') as fh:
    json.dump(out, fh)
print('field written:', out_path)
