@echo off
REM Drag and drop STL or 3MF files (or a folder of them) onto this file.
REM STEP files are written next to the originals.
REM Runs the copy next to this file, whether or not stl2curves is pip-installed.
set "PYTHONPATH=%~dp0;%PYTHONPATH%"
python -m stl2curves %*
echo.
pause
