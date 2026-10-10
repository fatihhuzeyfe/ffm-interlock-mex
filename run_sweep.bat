@echo off
REM ---------------------------------------------------------------------------
REM  Full FE sweep on the university PC. Nothing to install: Abaqus' own Python
REM  runs the launcher. Edit the two lines below, then double-click or run from
REM  a command prompt in this folder. Safe to stop (Ctrl+C) and start again:
REM  finished analyses are skipped.
REM ---------------------------------------------------------------------------
set ABQ=abaqus
set PAR=6

cd /d "%~dp0"
call %ABQ% python scripts\preflight.py --abaqus %ABQ%
echo.
echo === 1) Griffith benchmark (verification D2, 15 analyses)
call %ABQ% python scripts\launch_abaqus.py --only benchmark_griffith --parallel 1 --abaqus %ABQ%
echo === 2) main sweep, stage r0
call %ABQ% python scripts\launch_abaqus.py --only __r0 --parallel %PAR% --abaqus %ABQ%
echo === 3) main sweep, stage r90
call %ABQ% python scripts\launch_abaqus.py --only __r90 --parallel %PAR% --abaqus %ABQ%
echo === 4) pack the results for the laptop
call %ABQ% python scripts\transfer.py pack-results
call %ABQ% python scripts\launch_abaqus.py --status --abaqus %ABQ%
pause
