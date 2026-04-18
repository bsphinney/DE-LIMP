@echo off
:: ============================================================================
::  Launch_DE-LIMP_WSL.bat — One-click launcher for DE-LIMP on Windows via WSL2
:: ============================================================================
::
::  Runs DE-LIMP natively in a WSL2 Ubuntu distro. No Docker needed.
::  On first run, installs R + all packages inside WSL (20-30 min).
::  Subsequent runs start in ~30 seconds.
::
::  Why WSL over Docker:
::    - No .NET, no 9p msize tuning, no SSH key chmod gymnastics
::    - Bioconductor packages (limpa, ComplexHeatmap, MOFA2) compile cleanly
::    - Faster startup once installed
::
::  Requirements:
::    - Windows 10 2004+ or Windows 11
::    - WSL2 with Ubuntu (installed automatically if missing)
::
:: ============================================================================

setlocal enableextensions
title DE-LIMP (WSL)

echo.
echo  ============================================
echo    DE-LIMP Proteomics (WSL2 launcher)
echo  ============================================
echo.

:: ----------------------------------------------------------------------------
:: 1. Check WSL is installed
:: ----------------------------------------------------------------------------
where wsl >nul 2>&1
if errorlevel 1 (
    echo  ERROR: wsl.exe not found. You need Windows 10 2004+ or Windows 11.
    echo.
    echo  To install WSL2, open PowerShell as Administrator and run:
    echo      wsl --install
    echo  Then restart Windows and re-run this launcher.
    echo.
    pause
    exit /b 1
)

:: ----------------------------------------------------------------------------
:: 2. Check that an Ubuntu distro is installed
:: ----------------------------------------------------------------------------
wsl --list --quiet 2>nul | findstr /I "Ubuntu" >nul
if errorlevel 1 (
    echo  No Ubuntu distro found in WSL. Installing Ubuntu now...
    echo  ^(This opens a separate window — follow the prompts, then close it and re-run this launcher.^)
    echo.
    wsl --install -d Ubuntu
    echo.
    echo  Ubuntu install triggered. After it finishes setting up, re-run this launcher.
    pause
    exit /b 0
)

:: ----------------------------------------------------------------------------
:: 3. Copy the setup script into WSL (so it works even if this bat lives
::    on a Windows path WSL can't easily reach, like OneDrive)
:: ----------------------------------------------------------------------------
set "SETUP_SRC=%~dp0delimp_wsl_setup.sh"
if not exist "%SETUP_SRC%" (
    echo  ERROR: delimp_wsl_setup.sh not found next to this .bat file.
    echo  Expected at: %SETUP_SRC%
    echo.
    echo  Make sure you cloned the full DE-LIMP repo, not just this file.
    pause
    exit /b 1
)

echo  Copying setup script into WSL...
:: Convert Windows path to WSL path via wslpath, then cp to ~
wsl -d Ubuntu -e bash -c "cp \"$(wslpath -u '%SETUP_SRC%')\" ~/delimp_wsl_setup.sh && chmod +x ~/delimp_wsl_setup.sh"
if errorlevel 1 (
    echo  ERROR: Failed to copy setup script into WSL.
    pause
    exit /b 1
)

:: ----------------------------------------------------------------------------
:: 4. Run the app (setup is auto-triggered on first run)
:: ----------------------------------------------------------------------------
echo.
echo  Starting DE-LIMP...
echo  If this is your first run, expect ~20-30 min of R package installs.
echo  Subsequent runs take under 1 minute.
echo.
echo  When the terminal shows "Listening on http://0.0.0.0:3838",
echo  your browser will open automatically.
echo.

:: Open browser after a short delay (detached from the wsl process)
start "" /B cmd /c "timeout /t 90 /nobreak >nul && start http://localhost:3838"

:: Run inside WSL — this blocks until the app is stopped with Ctrl+C
wsl -d Ubuntu -e bash -c "bash ~/delimp_wsl_setup.sh"

echo.
echo  DE-LIMP has stopped.
pause
endlocal
