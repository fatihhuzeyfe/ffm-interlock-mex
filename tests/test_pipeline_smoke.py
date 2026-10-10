"""End-to-end smoke test of make_cases -> (synthetic FE data) -> run_postprocess on a real
interlock setup with the orthotropic literature materials. Checks plumbing, not mechanics:
every path is evaluated, onset loads are finite and positive, CSV files are written, and the
orientation-dependent bulk criterion is used."""
import sys, json, math, pathlib, subprocess, shutil
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _synthetic_raw(d):
    case = json.loads((d / "case.json").read_text())
    P = np.array(case["interface_points"], float)
    sing = np.array([s["xy"] for s in case["singular"]], float) if case.get("singular") else np.zeros((0, 2))
    seg_len = np.hypot(*np.diff(P, axis=0).T)
    S = np.concatenate([[0], np.cumsum(seg_len)])

    def field(x, y):                       # mild singular field around the corners
        r = np.min(np.hypot(sing[:, 0] - x, sing[:, 1] - y)) if len(sing) else 1.0
        return 1.0 + 0.4 * (r + 0.02) ** -0.3

    raw = d / "raw"
    raw.mkdir(exist_ok=True)
    rows = []
    s_nodes = np.unique(np.concatenate([np.linspace(0, S[-1], 800), [b[2] for b in case["interface_breakpoints"]]]))
    for s_ in s_nodes:
        seg = int(min(np.searchsorted(S, s_, side="right") - 1, len(P) - 2))
        f = (s_ - S[seg]) / max(seg_len[seg], 1e-12)
        x, y = P[seg] + f * (P[seg + 1] - P[seg])
        v = field(x, y)
        st = [0.2 * v, v, 0.1 * v]
        rows.append([float(s_), seg, float(x), float(y), st, st, [0.0, 0.01 * v]])
    bulk = {}
    for bl in case["bulk_lines"]:
        o = np.array(bl["origin"]); a = math.radians(bl["angle_deg"])
        ss = np.unique(np.concatenate([np.linspace(0, bl["L"], 200), [np.hypot(*(np.array(p) - o)) for p in bl["points"]]]))
        bulk[bl["id"]] = [[float(s_), 0, 0.3 * field(*(o + s_ * np.array([math.cos(a), math.sin(a)]))),
                           field(*(o + s_ * np.array([math.cos(a), math.sin(a)]))), 0.0, 0.0, 0.01] for s_ in ss]
    U0 = 0.07
    (raw / "uncracked.json").write_text(json.dumps(dict(ALLSE=U0, ALLWK=U0, F=8.0, interface=rows, bulk=bulk)))
    for p in case["paths"]:
        for k, cr in enumerate(p["cracks"]):
            A = cr["da"] if p["mode"] == "edge" else 2 * cr["da"]
            G = 2e-4 * (1 + 0.5 * math.sqrt(cr["da"]))
            if p["kind"] == "interface":
                lo, hi = cr["S_interval"]
                J = [[float(r[0]), r[1], 0.0, 2 * G * A / 0.01 / 20] for r in rows if lo <= r[0] <= hi][:20]
            else:
                J = [[float(r[0]), 0, 0.0, 2 * G * A / 0.01 / 20] for r in bulk[p["pid"]] if r[0] <= cr["da"]][:20]
            (raw / f"{p['pid']}__{k:02d}.json").write_text(json.dumps(dict(ALLSE=U0 - G * A, ALLWK=0, F=8.0, jumps=J)))


def test_pipeline_orthotropic(tmp_path):
    runs = tmp_path / "runs"
    study = json.loads((ROOT / "config" / "study.json").read_text())
    study["runs_dir"] = str(runs)
    study["post"] = dict(k_sigma=[1.0, 0.5], k_G=[1.0], feature_width_mm=[2.0], GIIc_over_GIc=[2.0])
    sp = tmp_path / "study.json"; sp.write_text(json.dumps(study))
    subprocess.run([sys.executable, str(ROOT / "scripts" / "make_cases.py"), "--study", str(sp),
                    "--only", "dovetail_ar1_a20__PLA-PETG__tension"], check=True, capture_output=True)
    dirs = sorted(runs.glob("dovetail_ar1_a20__PLA-PETG__tension__r*"))
    assert [d.name[-3:] for d in dirs] == ["_r0", "r90"]
    case = json.loads((dirs[0] / "case.json").read_text())
    assert "ortho" in case["fe"]["materials"]["A"] and len(case["paths"][0]["da"]) == study["da"]["n"]
    for d in dirs:
        _synthetic_raw(d)
    out = tmp_path / "res"
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "run_postprocess.py"), "--study", str(sp),
                        "--only", "dovetail_ar1_a20", "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    import csv
    rows = list(csv.DictReader(open(out / "paths.csv")))
    assert len(rows) == 2 * 2 * len(case["paths"])
    lam = [float(x["lam_c"]) for x in rows]
    assert all(np.isfinite(lam)) and min(lam) > 0
    assert (out / "governing.csv").exists() and (out / "ranking.csv").exists()


def test_bulk_orientation_dependence():
    from ffm.materials import bulk_props_for_crack, load_materials
    m = load_materials(ROOT / "config" / "materials.json")["materials"]["PLA"]
    along = bulk_props_for_crack(m, 0.0, 0.0)      # crack along the roads
    across = bulk_props_for_crack(m, 90.0, 0.0)    # crack across the roads
    assert math.isclose(along[0], m["Y"]) and math.isclose(across[0], m["X"])
    assert math.isclose(along[1], m["G_L"]) and math.isclose(across[1], m["G_T"])
    # rotating roads and crack together changes nothing
    assert np.allclose(bulk_props_for_crack(m, 120.0, 90.0), bulk_props_for_crack(m, 30.0, 0.0))
