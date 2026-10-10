@echo off
REM ---------------------------------------------------------------------------
REM  Cohesive-zone (CZM) runs: validation of the FFM onset + post-onset response.
REM      run_czm.bat          all cases (cases are generated if runs_czm does not exist)
REM      run_czm.bat test     one case only (check the model first)
REM ---------------------------------------------------------------------------
set PY=python
set ABQ=abaqus
set PAR=20
cd /d "%~dp0"
title ffm sweep czm
if not exist runs_czm %PY% scripts\make_czm_cases.py
if "%1"=="test" (
  %PY% scripts\launch_abaqus.py --runs runs_czm --only dovetail_ar1_a20__PLA-PETG__tension__r0__ks1__kg1__w5 --parallel 1 --abaqus %ABQ%
  goto :eof
)
%PY% scripts\launch_abaqus.py --runs runs_czm --parallel %PAR% --abaqus %ABQ%
%PY% scripts\czm_postprocess.py
pause
