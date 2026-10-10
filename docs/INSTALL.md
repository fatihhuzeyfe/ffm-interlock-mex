# Installation

## Python part (all platforms)

Python 3.9 to 3.12 is supported. From the repository root:

```bash
python -m venv .venv
# Windows (cmd / PowerShell):  .venv\Scripts\activate
# Linux / macOS:               source .venv/bin/activate
pip install -e ".[dev]"
python -m pytest -q
```

`pip install -e .` installs the `ffm` package in editable mode, so the scripts in `scripts/` and
`analysis/` can be run from any folder. The `[dev]` extra adds `pytest`.

With conda:

```bash
conda env create -f environment.yml
conda activate ffm-interlock
pip install -e . --no-deps
```

Dependencies: NumPy, SciPy, pandas, matplotlib, scikit-learn (surrogate models and Sobol
analysis), pytest.

## Abaqus (only for new finite element analyses)

* Abaqus/Standard 2017 or newer. The scripts in `abaqus/` are executed by Abaqus' own Python
  (Python 2.7 in 2017-2020, Python 3 in 2021 and later) and are written for both.
* The command that starts Abaqus must work in a terminal, e.g. `abaqus information=release`.
  If your installation uses another name (`abq2024`, `abaqus2017`), pass it with
  `--abaqus abq2024` to the launcher and set `ABQ=...` in the batch files.
* Check the installation and the license before a long run:

```bash
abaqus python scripts/preflight.py --abaqus abaqus
python scripts/launch_abaqus.py --check --parallel 1 --only flat__PLA-PETG__tension__r0
```

## MATLAB (only for the figures of the article)

MATLAB R2016b or newer, no toolboxes. The figure scripts read the CSV/JSON files written to
`figures/data/` in step 2 of [REPRODUCE.md](REPRODUCE.md) and write EPS, PDF and PNG files to
`figures/output_matlab/`.
