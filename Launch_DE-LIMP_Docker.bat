@echo off
title DE-LIMP Proteomics
color 1F

echo.
echo  ============================================
echo    DE-LIMP Proteomics - Starting...
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

:: Pull latest code
cd /d "%~dp0"
echo  Updating to latest version...
git pull --ff-only 2>nul
echo.

:: Auto-setup SSH key for HPC search (one-time)
if not exist "data\ssh" mkdir "data\ssh"
if not exist "data\ssh\id_ed25519" (
    :: Try common locations
    if exist "%USERPROFILE%\.ssh\id_ed25519" (
        echo  Found SSH key at %USERPROFILE%\.ssh\id_ed25519
        copy "%USERPROFILE%\.ssh\id_ed25519" "data\ssh\id_ed25519" >nul
        echo  Copied to data\ssh\ for HPC access.
        echo.
    ) else (
        echo  No SSH key found. To use HPC search features:
        echo    1. Place your id_ed25519 key in the data\ssh\ folder
        echo    2. Or generate one: ssh-keygen -t ed25519
        echo    3. Upload the public key to hippo.ucdavis.edu
        echo.
    )
)

:: Create data dirs if needed
if not exist "data\raw" mkdir "data\raw"
if not exist "data\fasta" mkdir "data\fasta"
if not exist "data\output" mkdir "data\output"

:: Build and start
echo  Starting DE-LIMP (this may take a minute)...
docker compose up -d --build >nul 2>&1

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
echo    URL:  http://localhost:3838
echo    SSH:  Key auto-loaded from data\ssh\
echo.
echo    For HPC search: click New Search tab,
echo    select Remote (SSH), click Test Connection.
echo.
start http://localhost:3838
echo  Press any key to stop DE-LIMP...
pause >nul

:: Cleanup
echo  Stopping DE-LIMP...
docker compose down >nul 2>&1
echo  Done!
timeout /t 2 >nul
