# =============================================================================
# build_diann_docker.ps1 — Build DIA-NN Docker image for use with DE-LIMP
# =============================================================================
#
# DIA-NN is proprietary software by Vadim Demichev. It is free for academic
# use but CANNOT be redistributed. This script downloads DIA-NN directly from
# the official GitHub release and builds a local Docker image.
#
# License: https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md
# Citation: Demichev V et al. (2020) Nature Methods 17:41-44
#
# Usage (PowerShell):
#   .\build_diann_docker.ps1              # Builds diann:2.0 (default)
#   .\build_diann_docker.ps1 -Version 2.1.0
#
# =============================================================================

param(
    [string]$Version = "2.0"
)

$ErrorActionPreference = "Stop"
$ImageName = "diann:$Version"
$DiannUrl = "https://github.com/vdemichev/DiaNN/releases/download/$Version/DIA-NN-$Version-Academia-Linux.zip"

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "  DIA-NN Docker Image Builder for DE-LIMP"
Write-Host "============================================================"
Write-Host ""
Write-Host "  Version:  $Version"
Write-Host "  Image:    $ImageName"
Write-Host ""
Write-Host "  IMPORTANT: DIA-NN is free for academic use." -ForegroundColor Yellow
Write-Host "  By proceeding, you agree to the DIA-NN license terms at:"
Write-Host "  https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md"
Write-Host ""
Write-Host "  Citation: Demichev V, Messner CB, Vernardis SI, Lilley KS,"
Write-Host "  Ralser M. DIA-NN: neural networks and interference correction"
Write-Host "  enable deep proteome coverage in high throughput."
Write-Host "  Nature Methods. 2020;17(1):41-44."
Write-Host "============================================================"
Write-Host ""

# Check for Docker
try {
    docker version | Out-Null
} catch {
    Write-Host "ERROR: Docker is not installed or not running." -ForegroundColor Red
    Write-Host "Install Docker Desktop from: https://www.docker.com/products/docker-desktop/"
    Write-Host "Make sure Docker Desktop is running before retrying."
    exit 1
}

# Check Docker daemon
try {
    docker info | Out-Null
} catch {
    Write-Host "ERROR: Docker daemon is not running." -ForegroundColor Red
    Write-Host "Start Docker Desktop and try again."
    exit 1
}

# Create temp build directory
$BuildDir = Join-Path $env:TEMP "diann-docker-build-$(Get-Random)"
New-Item -ItemType Directory -Path $BuildDir -Force | Out-Null

try {
    Set-Location $BuildDir

    # Download
    $ZipFile = "DIA-NN-$Version-Academia-Linux.zip"
    Write-Host "Downloading DIA-NN $Version Linux release..."
    Write-Host "URL: $DiannUrl"
    Write-Host ""

    try {
        Invoke-WebRequest -Uri $DiannUrl -OutFile $ZipFile -UseBasicParsing
    } catch {
        Write-Host "ERROR: Download failed." -ForegroundColor Red
        Write-Host "Check that version $Version exists at:"
        Write-Host "  https://github.com/vdemichev/DiaNN/releases"
        exit 1
    }

    if (-not (Test-Path $ZipFile) -or (Get-Item $ZipFile).Length -eq 0) {
        Write-Host "ERROR: Downloaded file is empty or missing." -ForegroundColor Red
        exit 1
    }

    Write-Host "Download complete: $([math]::Round((Get-Item $ZipFile).Length / 1MB, 1)) MB"

    # Extract
    Write-Host "Extracting..."
    Expand-Archive -Path $ZipFile -DestinationPath "diann-extract" -Force

    # Find the diann-linux binary
    $DiannBin = Get-ChildItem -Path "diann-extract" -Recurse -Filter "diann-linux" | Select-Object -First 1
    if (-not $DiannBin) {
        Write-Host "ERROR: Could not find diann-linux binary in the archive." -ForegroundColor Red
        Write-Host "Archive contents:"
        Get-ChildItem -Path "diann-extract" -Recurse | ForEach-Object { Write-Host "  $_" }
        exit 1
    }

    $DiannDir = $DiannBin.DirectoryName
    Write-Host "Found DIA-NN binary at: $($DiannBin.FullName)"

    # Copy binaries to build context
    New-Item -ItemType Directory -Path "diann-bin" -Force | Out-Null
    Copy-Item -Path "$DiannDir\*" -Destination "diann-bin\" -Force

    # Write Dockerfile
    @"
FROM --platform=linux/amd64 debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    libgomp1 \
    libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q https://dot.net/v1/dotnet-install.sh -O /tmp/dotnet-install.sh && \
    chmod +x /tmp/dotnet-install.sh && \
    /tmp/dotnet-install.sh --channel 8.0 --install-dir /usr/share/dotnet && \
    ln -s /usr/share/dotnet/dotnet /usr/bin/dotnet && \
    rm /tmp/dotnet-install.sh

ENV DOTNET_ROOT=/usr/share/dotnet
ENV PATH="`$PATH:/usr/share/dotnet"

COPY diann-bin/ /opt/diann/

RUN chmod +x /opt/diann/diann-linux 2>/dev/null || true && \
    ln -sf /opt/diann/diann-linux /usr/local/bin/diann

ENV LD_LIBRARY_PATH="/opt/diann:`${LD_LIBRARY_PATH}"

RUN diann --help 2>&1 | head -3 || echo "Note: diann --help returned non-zero (may be normal)"

WORKDIR /work

ENTRYPOINT ["diann"]
"@ | Set-Content -Path "Dockerfile" -Encoding UTF8

    # Build
    Write-Host ""
    Write-Host "Building Docker image: $ImageName" -ForegroundColor Cyan
    Write-Host "Platform: linux/amd64 (x86_64)"
    Write-Host ""

    docker build --platform linux/amd64 -t $ImageName .

    if ($LASTEXITCODE -ne 0) {
        Write-Host "ERROR: Docker build failed." -ForegroundColor Red
        exit 1
    }

    Write-Host ""
    Write-Host "============================================================" -ForegroundColor Green
    Write-Host "  DIA-NN $Version Docker image built successfully!" -ForegroundColor Green
    Write-Host "============================================================" -ForegroundColor Green
    Write-Host ""
    Write-Host "  Image name: $ImageName"
    Write-Host ""
    Write-Host "  Test with:"
    Write-Host "    docker run --rm $ImageName --help"
    Write-Host ""
    Write-Host "  To use with DE-LIMP:"
    Write-Host "    1. Start DE-LIMP"
    Write-Host "    2. Go to 'New Search' tab"
    Write-Host "    3. Select 'Local (Docker)' backend"
    Write-Host "    4. The image '$ImageName' will be detected automatically"
    Write-Host ""

    # Quick test
    Write-Host "Running quick verification..."
    docker run --rm --platform linux/amd64 $ImageName --help 2>&1 | Out-Null
    if ($LASTEXITCODE -eq 0) {
        Write-Host "  DIA-NN responds to --help" -ForegroundColor Green
    } else {
        Write-Host "  Warning: --help returned non-zero (may be normal for some versions)" -ForegroundColor Yellow
    }

    Write-Host ""
    Write-Host "Done!" -ForegroundColor Green

} finally {
    # Cleanup
    Set-Location $env:USERPROFILE
    Remove-Item -Path $BuildDir -Recurse -Force -ErrorAction SilentlyContinue
}
