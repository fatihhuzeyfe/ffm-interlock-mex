"""Pre-flight check on a new computer. Standard library only, so it runs with Abaqus' own
Python (no installation needed):

    abaqus python scripts/preflight.py --abaqus abaqus

Checks: Python version, project path, Abaqus command and release, CPU cores (-> suggested
--parallel), free disk space, state of the runs folder.
"""
import os, sys, json, shutil, subprocess, argparse, pathlib, re

ROOT = pathlib.Path(__file__).resolve().parents[1]


def sh(cmd, timeout=120):
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout,
                           errors="ignore", cwd=str(ROOT))
        return r.returncode, (r.stdout or "") + (r.stderr or "")
    except Exception as e:                       # noqa
        return 1, str(e)


def physical_cores():
    if os.name == "nt":
        rc, out = sh('powershell -NoProfile -Command "(Get-CimInstance Win32_Processor | '
                     'Measure-Object -Property NumberOfCores -Sum).Sum"', 60)
        m = re.search(r"\d+", out or "")
        if rc == 0 and m:
            return int(m.group())
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abaqus", default="abaqus")
    a = ap.parse_args()
    ok = True
    print(f"Python {sys.version.split()[0]} ({sys.executable})")
    p = str(ROOT)
    print(f"project: {p}")
    if " " in p or not p.isascii():
        print("  WARNING: path contains spaces or non-ASCII characters (e.g. Turkish letters); "
              "Abaqus is safer in a path like C:\\work\\ffm_mex")
        ok = False
    rc, out = sh(f"{a.abaqus} information=release", 180)
    rel = re.search(r"Abaqus\s+(\d{4})", out or "")
    if rc != 0 or not rel:
        print(f"  ERROR: '{a.abaqus} information=release' failed. Try --abaqus abq2025 (or the full path "
              "to abaqus.bat in the Commands folder).")
        print("  output:", (out or "").strip()[:400])
        ok = False
    else:
        print(f"Abaqus release: {rel.group(1)}  (2024+ -> Python 3 inside Abaqus; the driver supports both)")
    logical = os.cpu_count() or 1
    phys = physical_cores()
    print(f"CPU: {logical} logical" + (f", {phys} physical cores" if phys else ""))
    base = phys or max(1, logical // 2)
    sys.path.insert(0, str(ROOT / "scripts"))
    from sysinfo import memory_gb
    tot, free = memory_gb()
    by_mem = 99
    if tot:
        by_mem = max(1, int((free - 4) // 2))      # ~2 GB per CAE + Standard process, 4 GB reserve
        print(f"RAM: {free:.1f} GB free of {tot:.1f} GB")
    print(f"  suggested --parallel {max(1, min(base - 1, by_mem, 16))} "
          "(cores - 1, about 2 GB RAM per process; reduce if the licence server limits tokens)")
    free = shutil.disk_usage(p).free / 1e9
    print(f"free disk: {free:.0f} GB" + ("  (WARNING: < 20 GB)" if free < 20 else ""))
    runs = ROOT / "runs"
    if runs.is_dir():
        n_set = n_done = n_exp = 0
        for c in runs.glob("*/case.json"):
            case = json.loads(c.read_text())
            n_set += 1
            n_exp += sum(len(q["cracks"]) for q in case["paths"]) + 1
            raw = c.parent / "raw"
            n_done += len(list(raw.glob("*.json"))) if raw.is_dir() else 0
        print(f"runs: {n_set} setups, {n_done}/{n_exp} FE analyses present")
    else:
        print("runs: none yet (generate cases on a machine with numpy: python scripts/make_cases.py, "
              "or copy the runs folder)")
        ok = False
    print("\nREADY" if ok else "\nfix the items above first")


if __name__ == "__main__":
    main()
