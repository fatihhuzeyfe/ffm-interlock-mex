# Onset predictions of the study

`governing.csv.gz`: all 67 500 onset predictions (300 FE setups x 225 combinations of interface
strength, toughness, feature width and mode ratio): predicted onset path, mechanism, onset stress,
crack length and quality flags. Columns are described in `docs/DATA.md`. Read with
`pandas.read_csv("data/governing.csv.gz")`. Regenerate with `scripts/run_postprocess.py`.
