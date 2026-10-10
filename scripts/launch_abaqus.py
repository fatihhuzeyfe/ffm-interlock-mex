"""Run the Abaqus driver for many setups in parallel (Windows or Linux).

    python scripts/launch_abaqus.py --only benchmark_griffith --check     # 1 model, writes .cae, stops
    python scripts/launch_abaqus.py --only dovetail_ar1_a20__PLA-PETG__tension --parallel 4
    python scripts/launch_abaqus.py --parallel 8                          # whole study (restartable)
    python scripts/launch_abaqus.py --status                              # progress only

Each process is  abaqus cae noGUI=abaqus/build_and_run.py -- <case.json> <slice> <n_slices>
A setup is split into --slices slices (default = --parallel for a single setup, 1 otherwise).
Finished raw/*.json files are skipped, so interrupting and relaunching is safe.
Every process uses Abaqus license tokens: choose --parallel accordingly.
"""
import sys, os, json, argparse, pathlib, subprocess, time
from concurrent.futures import ThreadPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from sysinfo import memory_gb, abaqus_processes, abaqus_process_list
DRIVER = ROOT / "abaqus" / "build_and_run.py"


def setups(runs, only):
    return sorted(p.parent for p in runs.glob("*/case.json") if only in p.parent.name)


def progress(d):
    case = json.loads((d / "case.json").read_text())
    n_exp = sum(len(p["cracks"]) for p in case["paths"]) + 1
    raw = d / "raw"
    n_done = 0
    if raw.is_dir():
        for f in raw.glob("*.json"):
            try:
                json.loads(f.read_text())
                n_done += 1
            except ValueError:
                pass                      # truncated file: the driver recomputes it
    return n_done, n_exp


def run_one(cmd, logfile, shell):
    t0 = time.time()
    with open(logfile, "a") as fh:
        fh.write("\n$ " + " ".join(cmd) + "\n")
        fh.flush()
        rc = subprocess.call(subprocess.list2cmdline(cmd) if shell else cmd, stdout=fh, stderr=subprocess.STDOUT,
                             cwd=str(ROOT), shell=shell)
    return rc, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="substring filter on setup names")
    ap.add_argument("--parallel", type=int, default=2, help="concurrent Abaqus processes")
    ap.add_argument("--slices", type=int, default=0, help="slices per setup (0 = auto)")
    ap.add_argument("--abaqus", default="abaqus", help="Abaqus command (e.g. abq2024)")
    ap.add_argument("--runs", default=str(ROOT / "runs"))
    ap.add_argument("--check", action="store_true",
                    help="uncracked model only + save .cae for visual inspection")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="start even if Abaqus processes are already running")
    ap.add_argument("--min-free-gb", type=float, default=4.0,
                    help="do not start a new Abaqus process while free RAM is below this")
    ap.add_argument("--stagger", type=float, default=8.0,
                    help="seconds between the first parallel start-ups (avoids CAE start-up races)")
    a = ap.parse_args()

    ds = [d for d in setups(pathlib.Path(a.runs), a.only)]
    if not ds:
        sys.exit("no setups found - run scripts/make_cases.py first")

    tot_d = tot_e = 0
    todo = []
    for d in ds:
        n_done, n_exp = progress(d)
        tot_d += n_done; tot_e += n_exp
        if n_done < n_exp or a.check:
            todo.append(d)
    print(f"{len(ds)} setups, {tot_d}/{tot_e} FE analyses done, {len(todo)} setups to run")
    if a.status:
        for d in ds:
            n_done, n_exp = progress(d)
            if n_done < n_exp:
                print(f"   {d.name:50s} {n_done:4d}/{n_exp}")
        return

    shell = os.name == "nt"          # abaqus.bat needs the shell on Windows
    n_running = abaqus_processes(mine_only=True)      # other users' Abaqus jobs are not ours
    if n_running and not a.force:
        names = ", ".join(f"{a} (pid {b})" for a, b, _ in abaqus_process_list(mine_only=True)[:8])
        sys.exit(f"{n_running} Abaqus processes are already running: {names}\n"
                 "Probably left over from a launcher that was closed (two drivers on the same setup collide).\n"
                 "Run stop_all.bat first, or add --force if they belong to someone else / are unrelated.")
    if a.slices:
        n_sl = a.slices
    else:
        n_sl = a.parallel if len(todo) == 1 else 1
    if a.check:
        n_sl = 1
    tasks = []
    for d in todo:
        for s in range(n_sl):
            cmd = [a.abaqus, "cae", f"noGUI={DRIVER}", "--", str(d / "case.json"), str(s), str(n_sl)]
            if a.check:
                cmd += ["uncracked_only", "save_cae"]
            tasks.append((d, s, cmd))
    print(f"launching {len(tasks)} Abaqus processes, {a.parallel} at a time")

    def work(it):
        i, t = it
        if i < a.parallel:
            time.sleep(i * a.stagger)
        waited = 0
        while True:                       # memory guard (no Task Manager needed)
            _, free = memory_gb()
            if free is None or free >= a.min_free_gb or waited > 3600:
                break
            time.sleep(30); waited += 30
        d, s, cmd = t
        rc, dt = run_one(cmd, d / f"launcher_s{s}.log", shell)
        n_done, n_exp = progress(d)
        tot, free = memory_gb()
        mem = f", RAM free {free:.0f}/{tot:.0f} GB" if tot else ""
        # with several slices per setup the others may still be running: judge only whole setups
        if n_done < n_exp and not a.check and n_sl == 1:
            mem += "   <-- INCOMPLETE: see abaqus_log.txt / launcher_s*.log"
        print(f"[{time.strftime('%H:%M:%S')}] {d.name} slice {s}: rc={rc}, {dt / 60:.1f} min, "
              f"{n_done}/{n_exp} done{mem}", flush=True)
        return rc

    with ThreadPoolExecutor(a.parallel) as ex:
        rcs = list(ex.map(work, list(enumerate(tasks))))
    bad = sum(1 for r in rcs if r != 0)
    incomplete = [d.name for d in todo if progress(d)[0] < progress(d)[1]]
    if incomplete and not a.check:
        print(f"{len(incomplete)} setups are incomplete (rerun the same command to retry them):")
        for n in incomplete[:30]:
            print("   ", n)
    print("finished;", "all processes returned 0" if not bad else f"{bad} processes returned errors "
          "(see runs/<setup>/abaqus_log.txt and launcher_s*.log)")
    if a.check and not bad:
        print("open runs/<setup>/uncracked_check.cae in Abaqus/CAE and inspect partitions, seams, mesh.")


if __name__ == "__main__":
    main()
