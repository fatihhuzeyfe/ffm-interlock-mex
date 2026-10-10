"""Completeness and integrity check of the FE results before moving to the next stage.

    python scripts/verify_runs.py --only __r0                 # report
    python scripts/verify_runs.py --only __r0 --fix           # + delete bad results / stale locks
    python scripts/verify_runs.py --runs runs_fine            # other runs folder

Per setup it checks
  * completeness: every expected result file present and valid JSON (truncated files)
  * uncracked model: positive strain energy and reaction force, interface and bulk data present
  * every cracked model: strain energy did not increase (displacement control), crack faces
    actually opened (seam assigned), energy/closure consistency of each path (0.95-1.05)
  * abaqus_log.txt: FAILED analyses (last error shown), stale .lck files in work folders
--fix deletes the result files that failed a check and stale lock files, so that simply
re-running the launcher recomputes exactly those analyses. The report is also written to
docs/verify_<runs>_<filter>.txt (send this one file instead of many logs).
"""
import sys, json, argparse, pathlib, re, time
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ffm.postprocess import Setup


def load(f):
    try:
        return json.loads(f.read_text())
    except (ValueError, OSError):
        return None


def check_setup(d, fix):
    issues, bad_files = [], []
    case = json.loads((d / "case.json").read_text())
    raw = d / "raw"
    expected = ["uncracked"] + [f"{p['pid']}__{k:02d}" for p in case["paths"] for k in range(len(p["cracks"]))]
    present = {f.stem: f for f in raw.glob("*.json")} if raw.is_dir() else {}
    missing = [e for e in expected if e not in present]
    for name, f in present.items():
        if load(f) is None:
            issues.append(f"truncated/invalid result file {f.name}")
            bad_files.append(f)
    unc = load(present["uncracked"]) if "uncracked" in present else None
    if unc is not None:
        if not unc.get("ALLSE", 0) > 0:
            issues.append("uncracked: strain energy <= 0"); bad_files.append(present["uncracked"])
        if abs(unc.get("F", 0)) <= 0:
            issues.append("uncracked: zero reaction force"); bad_files.append(present["uncracked"])
        if not unc.get("interface"):
            issues.append("uncracked: no interface data"); bad_files.append(present["uncracked"])
        ids = {b["id"] for b in case.get("bulk_lines", [])}
        if ids - set(unc.get("bulk", {})):
            issues.append(f"uncracked: bulk profiles missing for {sorted(ids - set(unc.get('bulk', {})))[:3]}")
            bad_files.append(present["uncracked"])
        U0 = unc.get("ALLSE", 0)
        for name, f in present.items():
            if name == "uncracked":
                continue
            r = load(f)
            if r is None:
                continue
            if r.get("ALLSE", 0) > U0 * (1 + 1e-5):
                issues.append(f"{name}: strain energy increased (crack not released?)"); bad_files.append(f)
            elif not r.get("jumps"):
                issues.append(f"{name}: no crack-face jumps (seam not assigned)"); bad_files.append(f)
        # energy / closure consistency per path
        try:
            st = Setup(d)
            for p in case["paths"]:
                pd = st.path_data(p)
                if pd is None:
                    continue
                cc = pd.cc_ratio[np.isfinite(pd.cc_ratio)]
                if len(cc) and not (0.95 < float(np.median(cc)) < 1.05):
                    issues.append(f"path {p['pid']}: closure/energy ratio {np.median(cc):.3f}")
        except Exception as e:                       # noqa
            issues.append(f"post-processing check failed: {e!r}"[:200])
    log = d / "abaqus_log.txt"
    if log.exists():
        txt = log.read_text(errors="ignore")
        # FAILED analyses that are still missing (later successful reruns are not reported)
        failed = [m for m in re.findall(r"(\S+) FAILED", txt) if m in missing]
        if failed:
            # report the error MESSAGE (solver ***ERROR lines / exception text), not the traceback
            msgs = []
            for name in sorted(set(failed)):
                i = txt.rfind(name + " FAILED")
                block = txt[i:i + 3000].split("\n[")[0]          # up to the next log entry
                lines = [l.strip() for l in block.splitlines()
                         if "***ERROR" in l or "Error:" in l or "error" in l.lower() and "File " not in l]
                msgs.append(f"{name}: " + (" | ".join(lines[:3])[:300] if lines else block[-200:].replace("\n", " ")))
            issues.append(f"{len(set(failed))} still-missing analyses FAILED:")
            issues.extend("      " + m for m in msgs[:6])
        if missing and "uncracked" in missing:
            tail = " | ".join(txt.strip().splitlines()[-4:])
            issues.append(f"log tail: {tail[:300]}")
    for lck in d.glob("work_s*/*.lck"):
        issues.append(f"stale lock file {lck.relative_to(d)}")
        bad_files.append(lck)
    if fix:
        for f in set(bad_files):
            try:
                f.unlink()
            except OSError:
                pass
    return len(expected), len(expected) - len(missing), missing, issues, len(set(bad_files))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--fix", action="store_true")
    a = ap.parse_args()
    runs = ROOT / a.runs
    dirs = sorted(p.parent for p in runs.glob("*/case.json") if a.only in p.parent.name)
    out = []
    tot_e = tot_d = n_ok = n_bad = 0
    to_rerun = []
    for d in dirs:
        n_exp, n_done, missing, issues, nb = check_setup(d, a.fix)
        tot_e += n_exp; tot_d += n_done; n_bad += nb
        if not missing and not issues:
            n_ok += 1
            continue
        to_rerun.append(d.name)
        out.append(f"\n{d.name}: {n_done}/{n_exp} present" + (f", {len(missing)} missing" if missing else ""))
        if missing:
            out.append(f"   missing: {', '.join(missing[:8])}" + (" ..." if len(missing) > 8 else ""))
        for i in issues[:12]:
            out.append("   " + i)
    head = [f"verify {a.runs} '{a.only}'  {time.strftime('%Y-%m-%d %H:%M')}",
            f"{len(dirs)} setups: {n_ok} complete and clean, {len(to_rerun)} need attention",
            f"FE analyses present: {tot_d}/{tot_e}"]
    if a.fix:
        head.append(f"--fix: {n_bad} bad result/lock files deleted -> rerun the launcher to recompute them")
    if to_rerun:
        head.append("rerun:  python scripts\\launch_abaqus.py --runs " + a.runs +
                    (f" --only {a.only}" if a.only else "") + " --parallel 20")
    else:
        head.append("ALL CLEAN - ready for the next stage")
    rep = "\n".join(head + out)
    print(rep)
    (ROOT / "docs").mkdir(exist_ok=True)
    tag = re.sub(r"[^A-Za-z0-9]+", "_", f"{a.runs}_{a.only}").strip("_")
    (ROOT / "docs" / f"verify_{tag}.txt").write_text(rep, encoding="utf-8")


if __name__ == "__main__":
    main()
