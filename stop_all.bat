@echo off
REM Stops the sweep cleanly: the launcher window AND every Abaqus process it started.
REM (Closing the launcher window or Ctrl+C alone leaves Abaqus processes running.)
REM Afterwards restart with:  run_all.bat sweep    - finished analyses are kept.
echo Stopping all Abaqus processes of this user ...
taskkill /F /T /FI "WINDOWTITLE eq ffm sweep*" >nul 2>&1
for %%P in (ABQcaeK.exe standard.exe pre.exe package.exe explicit.exe ABQSMAStdDriver.exe) do (
  taskkill /F /IM %%P >nul 2>&1
)
timeout /t 5 >nul
python scripts\sysinfo.py
echo Abaqus processes should now be 0. Restart with: run_all.bat sweep
pause
