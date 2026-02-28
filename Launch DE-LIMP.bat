@echo off
title DE-LIMP Launcher
set "PS1=%~dp0launch_delimp.ps1"
set "URL=https://raw.githubusercontent.com/bsphinney/DE-LIMP/main/launch_delimp.ps1"

if not exist "%PS1%" (
    echo Downloading launcher script from GitHub...
    PowerShell -Command "Invoke-WebRequest -Uri '%URL%' -OutFile '%PS1%' -UseBasicParsing"
)

PowerShell -ExecutionPolicy Bypass -File "%PS1%"
pause
