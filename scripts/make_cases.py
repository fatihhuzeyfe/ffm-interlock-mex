"""Generate runs/<setup>/case.json for every FE setup of the study (and the benchmark).

    python scripts/make_cases.py                   # full study from config/study.json
    python scripts/make_cases.py --only dovetail   # subset (substring match on setup name)
    python scripts/make_cases.py --benchmark       # Griffith FE verification case
"""
import sys, json, argparse, pathlib, itertools, math
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from ffm.geometry import (GeomSpec, interface_polyline, candidate_paths, mirror_prune, build_case,
                          default_da, arc_length, sine_centres)
from ffm.materials import load_materials, dundurs, ortho_fe_entry, ortho_check


def expand_shapes(shapes):
    for s in shapes:
        keys = [k for k in s if k != "shape"]
        vals = [v if isinstance(v, list) else [v] for v in (s[k] for k in keys)]
        for combo in itertools.product(*vals):
            yield dict(shape=s["shape"], **dict(zip(keys, combo)))


def contact_wanted(study, shape, load):
    """Frictionless crack-face contact for the setups selected in study["crack_contact"]
    (where closing crack faces were found to overlap, e.g. curved interfaces in shear)."""
    rule = study.get("crack_contact")
    if not rule:
        return False
    return load in rule.get("loads", []) and shape in rule.get("shapes", [])


def fe_block(study, mats, A, B, load, H, raster=None):
    """Elastic data for the driver. With study["material_model"] == "orthotropic" and
    complete engineering constants, each material gets an orthotropic block and its road
    direction raster_deg (global, from the x axis); otherwise isotropic (E, nu)."""
    ortho = study.get("material_model", "isotropic") == "orthotropic"

    def one(name):
        d = dict(name=name, E=mats[name]["E"], nu=mats[name]["nu"])
        oc = ortho_fe_entry(mats[name]) if ortho else None
        if oc:
            ortho_check(mats[name])
            d.update(ortho=oc, raster_deg=float(raster or 0.0))
        return d
    return dict(materials=dict(A=one(A), B=one(B)),
                load=load, plane=study["plane"], delta=study["nominal_strain"] * 2 * H,
                mesh=study["mesh"], save_cae=False, keep_odb=False, ncpus=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default=str(ROOT / "config/study.json"))
    ap.add_argument("--materials", default=str(ROOT / "config/materials.json"))
    ap.add_argument("--only", default=None)
    ap.add_argument("--benchmark", action="store_true")
    a = ap.parse_args()
    study = json.loads(pathlib.Path(a.study).read_text(encoding="utf-8"))
    mats = load_materials(a.materials)["materials"]
    runs = ROOT / study["runs_dir"]
    da = default_da(study["da"]["n"], study["da"]["min"], study["da"]["max"])
    gd = study["geometry_defaults"]
    manifest = []

    if a.benchmark:
        write_benchmark(runs, study, mats, manifest)
    else:
        for spec, (A, B), load in itertools.product(list(expand_shapes(study["shapes"])),
                                                    study["pairs"], study["loads"]):
            g = GeomSpec(**{**gd, **spec})
            try:
                P = interface_polyline(g)
            except ValueError as e:
                print("skip:", e); continue
            centred = sine_centres(g, P) if g.shape == "sine" else None
            ext = None
            if study["da"].get("extend"):
                ext = dict(n=study["da"]["n"], min=study["da"]["min"], max=study["da"]["max_ext"],
                           max_bulk=study["da"].get("max_ext_bulk", study["da"]["max_ext"]))
            periodic = study.get("lateral", "free") == "periodic"
            paths = candidate_paths(g, P, da, bulk_offsets=study["bulk_offsets_deg"], centred_S=centred,
                                    extend=ext, periodic=periodic)
            if load == "tension" and study["symmetry_prune_for_tension"] and g.n_features == 1:
                paths = mirror_prune(paths, P, g.W)
            case = build_case(g, paths, P)
            case["lateral"] = study.get("lateral", "free")
            rasters = study.get("raster_deg", [None]) if study.get("material_model") == "orthotropic" else [None]
            for r in rasters:
                name = f"{g.tag()}__{A}-{B}__{load}" + ("" if r is None else f"__r{int(r)}")
                if a.only and a.only not in name:
                    continue
                c = dict(case)
                c.update(name=name, project_root=str(ROOT), fe=fe_block(study, mats, A, B, load, case["H"], r))
                if contact_wanted(study, g.shape, load):
                    c["fe"]["crack_contact"] = True
                write(runs, c, manifest)
    (runs / ("manifest_benchmark.json" if a.benchmark else "manifest.json")).write_text(json.dumps(manifest, indent=1))
    n = sum(m["n_cracked"] + 1 for m in manifest)
    print(f"{len(manifest)} FE setups, {n} FE analyses in total -> {runs}")


def _canon(case):
    # round-trip through JSON first: the in-memory case holds tuples and numpy floats, the file
    # holds lists and plain floats - without this the comparison reports false changes
    case = json.loads(json.dumps(case, default=float))

    def r(x):
        if isinstance(x, float):
            return float(f"{x:.9g}")
        if isinstance(x, list):
            return [r(v) for v in x]
        if isinstance(x, dict):
            return {k: r(v) for k, v in x.items() if k != "project_root"}
        return x
    return json.dumps(r(case), sort_keys=True)


def write(runs, case, manifest):
    d = runs / case["name"]
    d.mkdir(parents=True, exist_ok=True)
    f = d / "case.json"
    if f.exists() and (d / "raw").is_dir():
        # the FE definition changed (e.g. contact added): keep the old results aside so that
        # they are never mixed with the new model
        if _canon(json.loads(f.read_text())) != _canon(case):
            k = 1
            while (d / f"raw_superseded_{k}").exists():
                k += 1
            (d / "raw").rename(d / f"raw_superseded_{k}")
            print(f"   {case['name']}: definition changed -> old results moved to raw_superseded_{k}")
    f.write_text(json.dumps(case))
    nc = sum(len(p["cracks"]) for p in case["paths"])
    manifest.append(dict(name=case["name"], dir=str(d), n_paths=len(case["paths"]), n_cracked=nc))


def write_benchmark(runs, study, mats, manifest):
    """Centre-cracked homogeneous plate, symmetric FFM growth from a pre-crack of half length a=1.
    Reference solution: scripts/benchmark_griffith.py (exact FFM for the Griffith crack)."""
    a = 1.0
    g = GeomSpec(shape="flat", w=1.0, pitch=40.0, H=20.0, n_features=1)
    P = interface_polyline(g)
    Sc = g.W / 2
    da = list(np.geomspace(0.01, 2.0, 14))
    paths = candidate_paths(g, P, da, centred_S=[(Sc, a)], precrack=(Sc - a, Sc + a), L_max=2.0)
    case = build_case(g, paths, P, precrack=(Sc - a, Sc + a))
    m = "PLA"
    case.update(name="benchmark_griffith", project_root=str(ROOT), benchmark=dict(a=a))
    fe = fe_block(dict(study, material_model="isotropic"), mats, m, m, "tension", case["H"])   # exact solution is isotropic
    fe["mesh"] = dict(h_global=2.0, h_line=0.02, n_per_edge=6)
    fe["save_cae"] = True
    # sliding grips: uniform uniaxial stress far from the crack, so sigma_ref = F/W is the stress
    # the crack actually sees (clamped grips on a square plate raise the centre stress by ~8 %)
    fe["grips"] = "sliding"
    case["fe"] = fe
    write(runs, case, manifest)


if __name__ == "__main__":
    main()
