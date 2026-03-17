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

:: Create data dirs if needed
if not exist "data\raw" mkdir "data\raw"
if not exist "data\fasta" mkdir "data\fasta"
if not exist "data\output" mkdir "data\output"

:: Check for SSH key
set "SSH_DIR=%USERPROFILE%\.ssh"
if exist "%SSH_DIR%\id_ed25519" (
    echo  SSH key found: %SSH_DIR%\id_ed25519
    echo  HPC search will be available.
) else if exist "%SSH_DIR%\id_rsa" (
    echo  SSH key found: %SSH_DIR%\id_rsa
    echo  HPC search will be available.
) else (
    echo  No SSH key found in %SSH_DIR%
    echo  HPC search will not be available.
    echo  To set up: ssh-keygen -t ed25519
    echo.
)
echo.

:: Generate docker-compose override with current user's SSH dir
:: This lets each user on a shared PC use their own SSH key
echo services: > docker-compose.override.yml
echo   delimp: >> docker-compose.override.yml
echo     volumes: >> docker-compose.override.yml
echo       - ./data:/data >> docker-compose.override.yml
echo       - %SSH_DIR%:/home/shiny/.ssh:ro >> docker-compose.override.yml
echo     environment: >> docker-compose.override.yml
echo       - DELIMP_DATA_DIR=/data >> docker-compose.override.yml
echo       - DELIMP_SSH_KEY=/home/shiny/.ssh/id_ed25519 >> docker-compose.override.yml
echo       - DELIMP_SSH_USER=%USERNAME% >> docker-compose.override.yml

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
echo    URL:   http://localhost:3838
echo    User:  %USERNAME%
echo    SSH:   %SSH_DIR% mounted into container
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
del docker-compose.override.yml 2>nul
echo  Done!
timeout /t 2 >nul
