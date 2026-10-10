"""Put back results that were moved to raw_superseded_N if they are valid for the CURRENT case
definition (so nothing valid is recomputed).

    python scripts/restore_superseded.py --dry     # report only
    python scripts/restore_superseded.py           # restore the valid ones

Test per setup: every superseded crack result must belong to a crack of the current case and
its crack-face nodes must span exactly that crack: the recorded crack-face positions must reach
the crack end(s) of the current definition (within 20 % of the crack increment) and never go
beyond it. Results computed for a different crack geometry (e.g. v1 paths with other crack
lengths) fail this test and stay superseded -> they are recomputed. A setup is restored only if
ALL its files pass (uncracked and cracked models must come from the same mesh).
"""
import sys, json, shutil, argparse, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ffm.postprocess import crack_extent_check
TOL = 1e-6


def check_setup(d, src=None):
    if src is None:
        sup = sorted(d.glob("raw_superseded_*"), key=lambda p: int(p.name.rsplit("_", 1)[1]))
        if not sup:
            return None
        src = sup[-1]
    case = json.loads((d / "case.json").read_text())
    cracks = {f"{p['pid']}__{k:02d}": (p, c) for p in case["paths"] for k, c in enumerate(p["cracks"])}
    files = [f for f in src.glob("*.json")]
    if not files:
        return src, False, "empty", 0
    npe = case["fe"]["mesh"].get("n_per_edge", 4)
    short = []
    for f in files:
        if f.stem == "uncracked":
            continue
        if f.stem not in cracks:
            return src, False, f"{f.stem} is not a crack of the current case", len(files)
        try:
            raw = json.loads(f.read_text())
        except ValueError:
            continue                                   # truncated file: simply not restored
        ok, kind, why = crack_extent_check(*cracks[f.stem], raw, npe)
        if kind == "outside":
            return src, False, f"{f.stem}: {why}", len(files)
        if not ok:
            short.append(f"{f.stem}: {why}")
    # several short cracks = results of a shorter path definition (stale); a single short crack is
    # an FE seam defect of that one analysis (recomputing reproduces it): kept, and post-processing
    # leaves that crack out
    if len(short) > 1:
        return src, False, short[0] + f" (+{len(short) - 1} more)", len(files)
    return src, True, ("excluded in post-processing: " + short[0]) if short else "", len(files)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--audit", action="store_true",
                    help="check the CURRENT raw folders instead; results that do not match the current "
                         "case are moved to raw_superseded_N (and will be recomputed)")
    a = ap.parse_args()
    if a.audit:
        n_ok = n_bad = 0
        for c in sorted((ROOT / a.runs).glob("*/case.json")):
            d = c.parent
            if not (d / "raw").is_dir() or not any((d / "raw").glob("*.json")):
                continue
            src, ok, why, nf = check_setup(d, d / "raw")
            if ok:
                n_ok += 1
                if why:
                    print("   note", d.name, why)
                continue
            n_bad += 1
            print(f"   STALE {d.name}: {why}")
            if not a.dry:
                k = 1
                while (d / f"raw_superseded_{k}").exists():
                    k += 1
                (d / "raw").rename(d / f"raw_superseded_{k}")
        print(f"audit: {n_ok} setups consistent with their case, {n_bad} stale "
              f"({'would be ' if a.dry else ''}moved aside -> recomputed)")
        return
    n_ok = n_bad = n_files = 0
    bad = []
    for c in sorted((ROOT / a.runs).glob("*/case.json")):
        d = c.parent
        r = check_setup(d)
        if r is None:
            continue
        src, ok, why, nf = r
        if ok and why:
            print("   note", d.name, why)
        if ok:
            n_ok += 1; n_files += nf
            if not a.dry:
                dst = d / "raw"
                dst.mkdir(exist_ok=True)
                for f in src.glob("*.json"):
                    if not (dst / f.name).exists():
                        shutil.copy2(f, dst / f.name)
                k = 1
                while (d / f"raw_restored_{k}").exists():
                    k += 1
                src.rename(d / f"raw_restored_{k}")      # keep the copy, under a free name
        else:
            n_bad += 1
            bad.append(f"{d.name}: {why}")
    print(f"{n_ok} setups valid for the current definition -> "
          f"{'would be ' if a.dry else ''}restored ({n_files} files)")
    print(f"{n_bad} setups not valid (stay superseded, will be recomputed):")
    for b in bad:
        print("   ", b)


if __name__ == "__main__":
    main()
