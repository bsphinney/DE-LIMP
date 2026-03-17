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
echo    URL: http://localhost:3838
echo.
echo    To stop: close this window or run
echo    "docker compose down" in this folder
echo.
start http://localhost:3838
echo  Press any key to stop DE-LIMP...
pause >nul

:: Cleanup
echo  Stopping DE-LIMP...
docker compose down >nul 2>&1
echo  Done!
timeout /t 2 >nul
