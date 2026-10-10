# Reproducing the results of the article

The complete study is regenerated from the configuration in `config/` in three steps. The onset
predictions of the article are included in `data/governing.csv.gz` for comparison.

| Step | Input | Output | Needs | Time |
|---|---|---|---|---|
| 1 | `config/` | raw FE results of 51 132 analyses (`runs/`) | Abaqus/Standard | days |
| 2 | `runs/` | `results/*.csv`, `figures/data/` | Python | ~15 min (1 core) |
| 3 | `figures/data/` | Figures 1-5 and S1 | MATLAB | minutes |

Run all commands from the repository root with the virtual environment activated
(see [INSTALL.md](INSTALL.md)). Check the installation first with the commands of the README.

---

## Step 1: finite element analyses

Windows users can run the whole chain with the batch files (edit the settings at the top first):

```bat
run_sweep.bat                 :: FE sweep of all setups (restartable)
run_mesh_convergence.bat      :: three-level mesh study
run_czm.bat                   :: cohesive-zone counterpart
stop_all.bat                  :: stops all running Abaqus processes
```

The same steps on any platform:

```bash
python scripts/make_cases.py                          # runs/<setup>/case.json for all setups
python scripts/make_cases.py --benchmark              # Griffith benchmark case
abaqus python scripts/preflight.py --abaqus abaqus    # installation and license check
python scripts/launch_abaqus.py --check --parallel 1 --only flat__PLA-PETG__tension__r0
python scripts/launch_abaqus.py --parallel 8          # all setups; repeat after interruptions
python scripts/verify_runs.py                         # completeness, energy, closure, overlap
python scripts/check_fe_benchmark.py                  # FE Griffith benchmark
python scripts/mesh_convergence.py prepare --parallel 4
python scripts/make_czm_cases.py && python scripts/launch_abaqus.py --runs runs_czm --parallel 8
```

`--parallel` should not exceed the number of physical cores or the number of Abaqus license
tokens. The launcher keeps at least `--min-free-gb` (default 4 GB) of RAM free, staggers the
starts and skips every analysis whose result file already exists, so it can be stopped and
restarted at any time. In the study, 51 131 of 51 132 analyses completed; the missing one
(`U_ar0.5__PETG-ABS__shear__r0`, one crack length of path `I1b`) is reported in the
Supplementary Material.

---

## Step 2: onset predictions and figure data

```bash
python scripts/run_postprocess.py --jobs 8     # results/paths.csv, governing.csv, ranking.csv
python scripts/mesh_convergence.py analyze     # mesh convergence tables
python scripts/czm_postprocess.py              # cohesive-zone comparison
python analysis/export_figure_data.py          # figures/data/*.csv
```

Compare with the published predictions:

```bash
mkdir published && python -c "import gzip,shutil; shutil.copyfileobj(gzip.open('data/governing.csv.gz'), open('published/governing.csv','wb'))"
python scripts/compare_results.py published results
```

For the meshes and damage fields of Figures 5 and S1, run on the computer that holds the ODB files:

```bash
python scripts/export_figure_meshes.py         # figures/data/odb_*.json
```

---

## Step 3: figures

```matlab
cd figures/matlab
make_all_figures          % -> figures/output_matlab/Figure_1.eps/.pdf/.png ... Figure_S1_verification.*
```
