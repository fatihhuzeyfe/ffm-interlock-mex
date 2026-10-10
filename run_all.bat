@echo off
REM ===========================================================================
REM  Whole study on ONE workstation: cases -> FE (Abaqus) -> post-processing.
REM  Edit the settings, then run from a command prompt in this folder:
REM      run_all.bat            (= run_all.bat sweep: FE sweep, mesh convergence, post-processing)
REM      run_all.bat setup      (ONLY for a fresh project: tests + case generation + all stages)
REM      run_all.bat mesh       (start at a stage: setup | bench | sweep | mesh | post)
REM  Every stage is restartable: finished FE analyses are skipped.
REM ===========================================================================
set PY=python
set ABQ=abaqus
set PAR=6
REM optional backup folder for result zips (e.g. inside OneDrive); empty = no copy
set BACKUP=

cd /d "%~dp0"
title ffm sweep
if not "%1"=="" goto %1
REM default: never regenerate cases implicitly
goto sweep

:setup
echo ===== [setup] tests and case files
%PY% -m pytest -q tests
%PY% scripts\make_cases.py
%PY% scripts\make_cases.py --benchmark
%PY% scripts\mesh_convergence.py make
call %ABQ% python scripts\preflight.py --abaqus %ABQ%

:bench
echo ===== [bench] Griffith benchmark (verification D2)
%PY% scripts\launch_abaqus.py --only benchmark_griffith --parallel 1 --abaqus %ABQ%
if not exist docs mkdir docs
%PY% scripts\check_fe_benchmark.py > docs\benchmark_check.txt
type docs\benchmark_check.txt

:sweep
echo ===== [sweep] main sweep r0, then r90
%PY% scripts\launch_abaqus.py --only __r0 --parallel %PAR% --abaqus %ABQ%
call :backup
%PY% scripts\launch_abaqus.py --only __r90 --parallel %PAR% --abaqus %ABQ%
call :backup

:mesh
echo ===== [mesh] mesh convergence, coarse and fine levels (verification D3)
%PY% scripts\launch_abaqus.py --runs runs_mesh_coarse --parallel %PAR% --abaqus %ABQ%
%PY% scripts\launch_abaqus.py --runs runs_fine --parallel %PAR% --abaqus %ABQ%

:post
echo ===== [post] post-processing and paper outputs
%PY% scripts\run_postprocess.py --jobs %PAR% --plots
%PY% scripts\run_postprocess.py --study config\study_mesh_coarse.json --out results_mesh_coarse --jobs %PAR%
%PY% scripts\run_postprocess.py --study config\study_mesh_fine.json --out results_fine --jobs %PAR%
%PY% scripts\mesh_convergence.py analyze
call :backup
%PY% scripts\launch_abaqus.py --status --abaqus %ABQ%
echo ===== done. Results: results\  docs\mesh_convergence\  docs\benchmark_check.txt
pause
goto :eof

:backup
%PY% scripts\transfer.py pack-results
if not "%BACKUP%"=="" copy /y ..\ffm_results_*.zip "%BACKUP%" >nul
goto :eof
