@echo off
REM Double-click to open the stl2curves studio in your browser, or drop an STL file on this file.
REM Keep this window open while you use the page; close it to stop the studio.
REM Runs the copy next to this file, whether or not stl2curves is pip-installed.
set "PYTHONPATH=%~dp0;%PYTHONPATH%"
python -m stl2curves.studio %*
echo.
pause
