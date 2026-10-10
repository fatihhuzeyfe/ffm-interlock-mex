# Troubleshooting

If a problem is not listed here, open an issue with the command, the complete error message and
the output of `python scripts/sysinfo.py`. You can also give this file, [AGENTS.md](../AGENTS.md)
and the error message to an AI coding assistant.

## Python

| Symptom | Cause | Solution |
|---|---|---|
| `ModuleNotFoundError: No module named 'ffm'` | package not installed / wrong folder | `pip install -e .` in the repository root, activate the virtual environment |
| `No module named 'pandas'` / `'sklearn'` | dependencies missing | `pip install -e ".[dev]"` |
| `python` opens the Microsoft Store (Windows) | Python not on PATH | install Python from python.org with "Add to PATH", or use `py -3.11` |
| onset loads differ from `data/governing.csv.gz` | modified code or config | `git status`; restore `config/` and `ffm/`; rerun |
| `PermissionError` while writing results | file open in Excel or locked by OneDrive | close the file; keep the project outside synchronised folders |

## Abaqus

| Symptom | Cause | Solution |
|---|---|---|
| `'abaqus' is not recognized` | Abaqus command not on PATH | use the full command name (`abq2024`) with `--abaqus` / `ABQ=` |
| analyses wait or stop with license errors | not enough tokens | lower `--parallel`; check with `abaqus licensing ru` |
| contact definition error in Abaqus 2022+ / 2017 | master/slave vs main/secondary keywords | use the driver of this repository (detects the version); do not edit the contact block |
| `rc != 0` for one setup, others fine | element distortion or contact non-convergence in one crack | `python scripts/verify_runs.py --only <setup>` shows the missing analysis; rerun the launcher, it retries only missing files |
| PC becomes unresponsive | too many parallel processes | lower `--parallel`; the launcher keeps `--min-free-gb` free (default 4) |
| a run was interrupted | power cut, reboot | rerun the same command; finished analyses are skipped |
| stray Abaqus processes after Ctrl+C | child processes keep running | `stop_all.bat` (Windows) or `pkill -f abaqus` |
| very slow file access | antivirus or cloud sync scanning work files | exclude the project folder, keep it outside OneDrive/Dropbox |

## MATLAB figures

| Symptom | Cause | Solution |
|---|---|---|
| raw TeX text such as `\ell` in labels | command not supported by the MATLAB TeX interpreter | use the scripts of this repository (only supported commands are used) |
| fonts differ in EPS | Helvetica substituted | `figures/matlab/ffm_style.m` sets the font; Arial is equivalent |
| `File not found ... odb_*.json` | Abaqus exports not present | the panels are drawn without them; run `python scripts/export_figure_meshes.py` where the ODBs are |
| CSV columns look wrong after opening in Excel | decimal comma locale | do not save the CSV from Excel; the scripts read the original files |
