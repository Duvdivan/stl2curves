@echo off
REM Drag and drop STL files (or a folder of them) onto this file.
REM STEP files are written next to the originals.
python "%~dp0stl2curves.py" %*
echo.
pause
