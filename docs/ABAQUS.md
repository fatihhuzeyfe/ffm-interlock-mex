# Finite element model and Abaqus driver

## Model

* Periodic unit cell of width `pitch * w` (default 4w) and height `2 H w` (default 16w); the
  geometry is defined in units of the feature width w (w = 1 in the FE model).
* Plane strain, linear elastic, transversely isotropic materials with axis 1 along the roads;
  raster 0 deg = roads parallel to the mean interface, 90 deg = perpendicular.
* Periodic constraint equations between the lateral edges (reference node carries the
  macroscopic stretch and shear); lower edge clamped; displacement of the upper edge
  (nominal strain 0.001) in tension (u2) or shear (u1).
* Quadratic plane-strain triangles (CPE6); element size `mesh.h_line` along all candidate
  crack paths (at least `mesh.n_per_edge` elements per crack increment), `mesh.h_global` elsewhere.
* Every crack is a seam in an otherwise identical mesh; frictionless contact on the crack faces
  for the shapes and loads listed in `crack_contact`.

## Driver

`abaqus/build_and_run.py` is started by `scripts/launch_abaqus.py` as

```
abaqus cae noGUI=abaqus/build_and_run.py -- runs/<setup>/case.json <slice> <n_slices>
```

It builds the uncracked model and one model per crack, runs them, and writes one JSON file per
analysis to `runs/<setup>/raw/` (format in [DATA.md](DATA.md)). Work files go to
`runs/<setup>/work_s<slice>/` and are deleted after extraction unless `keep_odb` is set.

Version compatibility: the contact API changed in Abaqus 2022 (master/slave -> main/secondary);
the driver detects the version and uses the right keywords. Results of 2017 and 2025 differ by
less than 4e-7 in the onset stress.

## Study definition (`config/study.json`)

| Key | Meaning |
|---|---|
| `shapes` | geometry families and parameters: `flat`; `sine`, `U`, `T` with aspect ratios `ar` = d/w; `dovetail` with `ar` and flank angles `alpha` |
| `pairs` | material pairs [A, B] (A = lower material) |
| `loads` | `tension`, `shear` |
| `raster_deg` | road direction(s) of both materials |
| `geometry_defaults` | `pitch`, `H`, `head_ratio`, `neck_frac`, `n_features` |
| `mesh` | `h_global`, `h_line`, `n_per_edge` |
| `da` | crack lengths: `n` per path, `min`, `max`, extension to the next feature (`max_ext`, `max_ext_bulk`) |
| `bulk_offsets_deg` | directions of the bulk paths relative to the bisector of each re-entrant corner |
| `crack_contact` | loads and shapes with crack-face contact |
| `post` | interface quality factors `k_sigma`, `k_G`, feature widths in mm, `GIIc_over_GIc` |
| `_version` | the paper uses `v2-final (frozen)`; change the version string when you change the study |

Changing a key changes the case definitions; `make_cases.py` detects this and marks only the
affected analyses for re-running.

## Materials (`config/materials.json`)

`materials` holds the bulk data (E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, strengths X, Y, S,
toughnesses G_L along and G_T across the roads); `interfaces` holds the reference interface of
each pair (rule `weaker_intralayer`: strength and toughness of the weaker road-to-road bond,
`GIIc_over_GIc`, energy law). The values used in the paper are listed in [MATERIALS.md](MATERIALS.md).

## Useful tools

| Script | Purpose |
|---|---|
| `scripts/plot_geometries.py` | unit cells and candidate paths of all geometries |
| `scripts/launch_abaqus.py --status` | progress of all setups |
| `scripts/verify_runs.py` | completeness and consistency checks of the raw results |
| `scripts/restore_superseded.py --audit` | consistency of results with the current case definitions |
| `scripts/fe_figure.py` | mesh and stress field of one analysis (needs the kept ODB) |
| `scripts/export_figure_meshes.py` | meshes and CZM damage for the figures |
| `scripts/sysinfo.py` | Python, Abaqus, CPU and memory information for bug reports |
