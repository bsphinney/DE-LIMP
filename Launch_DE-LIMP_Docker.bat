@echo off
title DE-LIMP Proteomics
color 1F

echo.
echo  ============================================
echo    DE-LIMP Proteomics
echo  ============================================
echo.

:: Check if Docker is running
docker info >nul 2>&1
if errorlevel 1 (
    echo  Docker Desktop is not running!
    echo.
    echo  Please start Docker Desktop and wait for it
    echo  to finish loading, then try again.
    echo.
    pause
    exit /b 1
)

cd /d "%~dp0"

:: Ask which HIVE user
echo  Who is using DE-LIMP today?
echo.

:: List available SSH keys in data/ssh/
if not exist "data\ssh" mkdir "data\ssh"
set "found_keys=0"
for %%f in (data\ssh\*) do (
    if /i not "%%~nxf"==".gitkeep" (
        echo    %%~nxf
        set "found_keys=1"
    )
)
if "%found_keys%"=="0" (
    echo    No SSH keys found yet.
    echo.
    echo  To set up HPC access, each user should:
    echo    1. Copy their SSH private key to:
    echo       %CD%\data\ssh\
    echo    2. Rename it to their HIVE username
    echo       Example: %CD%\data\ssh\jsmith
    echo.
)
echo.
set /p HIVE_USER="  Enter your HIVE username (or press Enter to skip HPC): "

:: Set up SSH key for this user
set "SSH_KEY_PATH="
if defined HIVE_USER (
    if exist "data\ssh\%HIVE_USER%" (
        set "SSH_KEY_PATH=/data/ssh/%HIVE_USER%"
        echo  Using SSH key: data\ssh\%HIVE_USER%
    ) else if exist "data\ssh\id_ed25519" (
        set "SSH_KEY_PATH=/data/ssh/id_ed25519"
        echo  Using default SSH key: data\ssh\id_ed25519
    ) else (
        echo  No SSH key found for %HIVE_USER%.
        echo  Place your key at: %CD%\data\ssh\%HIVE_USER%
    )
)
echo.

:: Pull latest code
echo  Updating to latest version...
git pull --ff-only 2>nul

:: Create data dirs
if not exist "data\raw" mkdir "data\raw"
if not exist "data\fasta" mkdir "data\fasta"
if not exist "data\output" mkdir "data\output"

:: Pass user info as environment variables
set "DELIMP_SSH_USER=%HIVE_USER%"
set "DELIMP_SSH_KEY=%SSH_KEY_PATH%"

:: Start container (only rebuild if Dockerfile changed)
echo.
echo  Starting DE-LIMP...
docker compose up -d 2>nul
if errorlevel 1 (
    echo.
    echo  Building container - this may take a few minutes on first run...
    echo.
    docker compose up -d --build
)

:: Wait for app to respond
echo  Waiting for app to start...
:wait_loop
timeout /t 3 /nobreak >nul
curl -s -o nul -w "%%{http_code}" http://localhost:3838 2>nul | findstr "200 302" >nul
if errorlevel 1 goto wait_loop

:: Open browser
echo.
echo  ============================================
echo    DE-LIMP is running!
echo  ============================================
echo.
echo    URL:   http://localhost:3838
if defined HIVE_USER (
    echo    User:  %HIVE_USER%
)
echo.
start http://localhost:3838
echo  Press any key to stop DE-LIMP...
pause >nul

:: Cleanup
echo  Stopping DE-LIMP...
docker compose down >nul 2>&1
echo  Done!
timeout /t 2 >nul
