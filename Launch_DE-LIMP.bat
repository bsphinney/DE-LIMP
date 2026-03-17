@echo off
title DE-LIMP Launcher
set "PS1=%~dp0launch_delimp.ps1"
set "URL=https://raw.githubusercontent.com/bsphinney/DE-LIMP/main/launch_delimp.ps1"

echo Updating launcher from GitHub...
PowerShell -Command "Invoke-WebRequest -Uri '%URL%' -OutFile '%PS1%' -UseBasicParsing" 2>nul

PowerShell -ExecutionPolicy Bypass -File "%PS1%"
pause
