# Guide for AI coding assistants

This file gives an AI assistant (or a new contributor) the context needed to work on this
repository efficiently. Read it before changing code or diagnosing an error.

## Purpose

Prediction of crack onset at periodic interlocked interfaces between two 3D-printed polymers with
the coupled stress-energy criterion of finite fracture mechanics (FFM). Abaqus computes linear
elastic fields once; all material, interface-quality and size parameters are applied afterwards
in Python.

## Architecture

| Path | Role | Runs in |
|---|---|---|
| `ffm/kernel.py` | `Criterion`, `ffm_path`: stress/energy load factors and coupled minimum | CPython 3.9+ |
| `ffm/geometry.py` | interface polylines, singular points, candidate paths, periodic cell | CPython |
| `ffm/postprocess.py` | reads raw JSON, tractions along paths, G from energy and closure, `evaluate()` | CPython |
| `ffm/materials.py` | material data, orthotropic strength/toughness, Dundurs parameters, l_ch | CPython |
| `ffm/abq_utils.py` | helpers shared with the Abaqus driver (no NumPy) | Abaqus Python 2.7 and 3 |
| `abaqus/build_and_run.py` | builds, meshes, runs and extracts one setup (uncracked + every crack); CZM mode | Abaqus/CAE noGUI |
| `scripts/make_cases.py` | `config/*.json` -> `runs/<setup>/case.json` (idempotent, hash-based change detection) | CPython |
| `scripts/launch_abaqus.py` | parallel, restart-safe Abaqus launcher | CPython |
| `scripts/run_postprocess.py` | FFM for all setups and parameter combinations -> `results/*.csv` | CPython |
| `analysis/export_figure_data.py` | figure data -> `figures/data/` | CPython |
| `figures/matlab/` | final figures (MATLAB R2016b+, fixed cm layout, EPS) | MATLAB |

Data flow: `config` -> `runs/<setup>/case.json` -> Abaqus -> `runs/<setup>/raw/*.json` ->
`results/*.csv` -> `figures/data/*.csv` -> figures. Formats: `docs/DATA.md`.

## Conventions and invariants

* FE geometry in units of the feature width (w = 1); MPa; N/mm. The physical width enters only
  in post-processing (`post.feature_width_mm`).
* Load factor scaling: stresses ~ λ, energy release rates ~ λ². Never re-run Abaqus to change
  strengths, toughnesses, k factors, widths or G_IIc/G_Ic.
* `config/study.json` with `_version: "v2-final (frozen)"` reproduces the paper. Do not edit it in
  place; copy it and pass `--study`.
* Raw file names: `uncracked.json`, `<pid>__<kk>.json`. A file is "done" if it is valid JSON.
* Code in `abaqus/` and `ffm/abq_utils.py` must stay compatible with Python 2.7 (Abaqus 2017):
  no f-strings, no NumPy, no pathlib.
* Abaqus 2022+ uses main/secondary instead of master/slave in contact definitions; the driver
  handles both.
* All text in the repository is English.

## Commands

```bash
pip install -e ".[dev]"                               # install
python -m pytest -q                                   # tests (seconds)
python scripts/make_cases.py --benchmark && python scripts/benchmark_griffith.py   # FFM chain vs exact
python scripts/run_postprocess.py --only <substr> --plots
python scripts/verify_runs.py --only <substr>
python scripts/launch_abaqus.py --status
```

## Definition of done for a change

1. `python -m pytest -q` passes.
2. `python scripts/benchmark_griffith.py` prints `BENCHMARK PASSED`; changes that alter the results
   are checked against `data/governing.csv.gz` with `scripts/compare_results.py`.
3. Abaqus-side changes: one setup run with `--check` and one full setup, followed by
   `verify_runs.py`.
4. Documentation updated (`README.md`, `docs/`).

## Known pitfalls

* Long Abaqus runs on Windows: keep the project outside OneDrive; antivirus can lock work files.
* `--parallel` larger than the license tokens makes jobs queue, not fail.
* Periodic cell: a crack reaching the cell edge ends there (duplicated nodes tied to one partner);
  in CZM mode the seam continues across the edge (face-by-face pairing).
* Crack-face overlap without contact makes onset loads conservative; flagged in `governing.csv`.
* MATLAB TeX interpreter does not support `\ell`; use `\itl\rm_{ch}`.
