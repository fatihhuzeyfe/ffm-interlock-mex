# Verification

| Level | Command | Result in the article |
|---|---|---|
| Unit and pipeline tests | `python -m pytest -q` | all pass |
| Comparison with the study | `python scripts/compare_results.py published results` | onset predictions against `data/governing.csv.gz` |
| FFM kernel vs closed form | `python scripts/make_cases.py --benchmark` then `python scripts/benchmark_griffith.py` | Griffith crack, coupled criterion |
| FE chain vs closed form | `python scripts/check_fe_benchmark.py` | Ḡ within -2.3/+2.1 %, closure/energy 0.9998-1.0001, onset stress within 0.68 % |
| Run integrity | `python scripts/verify_runs.py` | 51 131 / 51 132 analyses; energy balance and closure checks |
| Mesh convergence | `python scripts/mesh_convergence.py analyze` | 9 setups x 225 combinations; median change 0.03 %, max 5.9 % (contact), GCI 0.12 % |
| Solver version | `python scripts/compare_results.py <results_2017> <results_2025>` | Abaqus 2017 vs 2025: < 4e-7 |
| CZM counterpart | `python scripts/czm_postprocess.py` | CZM/FFM 0.89-0.97 (flat), 0.77-0.95 (PLA/ABS), 0.50-0.92 (PLA/PETG) |

## What the checks do not show

* Agreement with experiments: the interface properties are scaled from a reference value and
  have to be measured for a given printer and process.
* Effects that are not modelled: three-dimensional effects, bulk fracture after debonding,
  friction, residual stresses, crack-face contact outside the listed cases.
* Known limitations of the data: one non-converged analysis; crack-face overlap in 7 % of the
  governing dovetail results under shear (conservative onset loads); flagged columns in
  `results/governing.csv` identify these cases.
