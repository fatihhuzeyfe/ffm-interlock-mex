@echo off
REM Coarse and fine meshes of the representative set (verification D3).
REM The case files must have been generated on the laptop with
REM   python scripts/mesh_convergence.py prepare   (then the make_cases lines it prints)
set ABQ=abaqus
set PAR=6
cd /d "%~dp0"
for %%L in (runs_mesh_coarse runs_fine) do (
  call %ABQ% python scripts\launch_abaqus.py --runs %%L --parallel %PAR% --abaqus %ABQ%
)
call %ABQ% python scripts\transfer.py pack-results
pause
