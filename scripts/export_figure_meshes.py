"""Export the Abaqus meshes and fields shown in the figures (run on the computer with the ODBs):

    python scripts/export_figure_meshes.py [--abaqus abaqus]

Writes figures/data/odb_*.json:
  odb_mesh_coarse / _medium / _fine   uncracked dovetail (d/w = 1, alpha = 20 deg, PLA/PETG, tension)
                                      at the three mesh levels, with the stress field
  odb_czm_<case>_t<time>              CZM runs: deformed mesh, stress and interface damage at
                                      selected step times (before / at / after onset)
"""
import sys, subprocess, argparse, pathlib, json

ROOT = pathlib.Path(__file__).resolve().parents[1]
SETUP = "dovetail_ar1_a20__PLA-PETG__tension__r0"
MESH = [("coarse", "runs_mesh_coarse"), ("medium", "runs"), ("fine", "runs_fine")]
CZM = [("dovetail_ar1_a20__PLA-PETG__tension__r0__ks1__kg1__w5", [0.12, 0.18, 0.30]),
       ("U_ar1__PLA-PETG__shear__r0__ks1__kg1__w5", [0.12, 0.25, 0.60]),
       ("flat__PLA-PETG__tension__r0__ks1__kg1__w5", [0.15, 0.25])]


def run(abq, odb, case, out, t=None):
    if not odb.exists():
        print("   missing:", odb); return
    cmd = f'{abq} python "{ROOT / "abaqus" / "export_figure_odb.py"}" "{odb}" "{case}" "{out}"'
    if t is not None:
        cmd += f" {t}"
    print("$", cmd)
    subprocess.call(cmd, shell=True, cwd=str(ROOT))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abaqus", default="abaqus")
    a = ap.parse_args()
    out = ROOT / "figures" / "data"
    out.mkdir(parents=True, exist_ok=True)
    for lvl, runs in MESH:
        d = ROOT / runs / SETUP
        run(a.abaqus, d / "work_s0" / "uncr.odb", d / "case.json", out / f"odb_mesh_{lvl}.json")
    for case, times in CZM:
        d = ROOT / "runs_czm" / case
        for t in times:
            run(a.abaqus, d / "work_s0" / "czm.odb", d / "case.json",
                out / f"odb_czm_{case.split('__')[0]}_{case.split('__')[2]}_t{int(round(100 * t)):03d}.json", t)
    print("done ->", out)


if __name__ == "__main__":
    main()
