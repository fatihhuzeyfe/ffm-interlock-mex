# Data formats

All lengths of the FE model are in units of the feature width w (w = 1); stresses in MPa;
energy release rates in N/mm (= kJ/m^2). Physical sizes enter only in post-processing.

## `runs/<setup>/case.json`

Complete definition of one FE setup, written by `scripts/make_cases.py`.

| Key | Content |
|---|---|
| `name`, `tag`, `geometry` | setup name (`<geometry>__<pair>__<load>__<raster>`), definition hash, geometry parameters |
| `W`, `H`, `lateral` | cell width, half height, lateral boundary condition (`periodic`) |
| `interface_points`, `interface_breakpoints` | interface polyline and its partition points |
| `singular` | singular points (corners) of the interface |
| `paths` | candidate crack paths: `pid`, `kind` (`interface`/`bulk`), polyline, `cracks` (crack lengths `da` and their position) |
| `bulk_lines` | straight bulk paths (origin, angle, length) |
| `fe` | materials (incl. raster), load, plane, nominal displacement `delta`, mesh |

Path identifiers: `I<n>f` / `I<n>b` = interface path from singular point n, forward/backward
along the interface; `C<n>` = centred interface crack (smooth interfaces); `B<n><mat><angle>` =
bulk path from corner n into material `A`/`B` at the given angle (degrees).

## `runs/<setup>/raw/*.json`

One file per analysis: `uncracked.json` and `<pid>__<k>.json` (crack length index k = 00...11).

| Key | Content |
|---|---|
| `ALLSE`, `ALLWK` | strain energy and external work of the model |
| `F` | reaction force of the loaded edge (per unit thickness) |
| `n_elem` | number of elements |
| `interface` (uncracked) | stresses along the interface: arc length, segment, x, y, stress tensor in A and B, displacement |
| `bulk` (uncracked) | stresses along every bulk path |
| `jumps` (cracked) | crack-face displacement jumps and nodal forces used for the crack-closure partition |

## `results/`

| File | Rows | Content |
|---|---|---|
| `paths.csv` | every path x 225 combinations | onset load of each candidate path (`lam_c`, `sigma_nom_c`, `da_c`, `governing` = stress/energy/coupled, flags) |
| `governing.csv` | 67 500 | weakest path of each setup and combination = predicted onset (`pid`, `kind`, `material`, `sigma_nom_c`, `da_c`, `mode_mix_at_c`, runner-up and margin, quality flags) |
| `ranking.csv` | 67 500 | `sigma_nom_c` relative to the flat interface (`rel_to_flat` = gain), onset location, rank |

Main columns: `k_sigma`, `k_G` (interface quality factors), `width_mm` (feature width w),
`GIIc_over_GIc`, `mu_interface` (= w / l_ch, the inverse of the interface length scale),
`sigma_nom_c` (onset nominal stress, MPa), `da_c` (crack length at onset, mm),
`segment_exhausted` (crack spans the whole interface segment), `interpenetration` /
`overlap_ratio` (crack-face overlap without contact; onset load is conservative).

## `figures/data/`

Not included in the repository; generated in step 2 of REPRODUCE.md. CSV files named after the figure and panel (`fig2a_master_curve_binned.csv`, ...), written by
`analysis/export_figure_data.py` and `scripts/czm_postprocess.py`; `odb_*.json` contain meshes,
stresses and interface damage exported from Abaqus by `scripts/export_figure_meshes.py`.

## `data/governing.csv.gz`

All 67 500 onset predictions of the study, same columns as `results/governing.csv` above.

## Raw results of the study

The raw results of all 51 132 analyses (about 0.5 GB) are not part of the repository; they are
regenerated with level 2 of [REPRODUCE.md](REPRODUCE.md) and are available from the author on
reasonable request.
