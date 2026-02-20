#!/bin/bash
# =============================================================================
# build_diann_docker.sh — Build DIA-NN Docker image for use with DE-LIMP
# =============================================================================
#
# DIA-NN is proprietary software by Vadim Demichev. It is free for academic
# use but CANNOT be redistributed. This script downloads DIA-NN directly from
# the official GitHub release and builds a local Docker image.
#
# License: https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md
# Citation: Demichev V et al. (2020) Nature Methods 17:41-44
#
# Usage:
#   bash build_diann_docker.sh [VERSION]
#
# Examples:
#   bash build_diann_docker.sh           # Builds diann:2.0 (default)
#   bash build_diann_docker.sh 2.0       # Same as above
#
# The resulting image can be used with DE-LIMP's "Local (Docker)" search mode.
# =============================================================================

set -euo pipefail

DIANN_VERSION="${1:-2.0}"
IMAGE_NAME="diann:${DIANN_VERSION}"

# DIA-NN Linux releases are distributed as zip files on GitHub
DIANN_URL="https://github.com/vdemichev/DiaNN/releases/download/${DIANN_VERSION}/diann-linux.zip"

echo "============================================================"
echo "  DIA-NN Docker Image Builder for DE-LIMP"
echo "============================================================"
echo ""
echo "  Version:  ${DIANN_VERSION}"
echo "  Image:    ${IMAGE_NAME}"
echo ""
echo "  IMPORTANT: DIA-NN is free for academic use."
echo "  By proceeding, you agree to the DIA-NN license terms at:"
echo "  https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md"
echo ""
echo "  Citation: Demichev V, Messner CB, Vernardis SI, Lilley KS,"
echo "  Ralser M. DIA-NN: neural networks and interference correction"
echo "  enable deep proteome coverage in high throughput."
echo "  Nature Methods. 2020;17(1):41-44."
echo "============================================================"
echo ""

# Check for Docker
if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker is not installed or not on PATH."
    echo "Install Docker Desktop from: https://www.docker.com/products/docker-desktop/"
    exit 1
fi

# Check Docker daemon is running
if ! docker info &> /dev/null; then
    echo "ERROR: Docker daemon is not running."
    echo "Start Docker Desktop and try again."
    exit 1
fi

# Create temporary build directory
BUILD_DIR=$(mktemp -d)
trap "rm -rf $BUILD_DIR" EXIT
cd "$BUILD_DIR"

echo "Downloading DIA-NN ${DIANN_VERSION} Linux release..."
if command -v wget &> /dev/null; then
    wget -q --show-progress "$DIANN_URL" -O diann-linux.zip
elif command -v curl &> /dev/null; then
    curl -L --progress-bar "$DIANN_URL" -o diann-linux.zip
else
    echo "ERROR: Neither wget nor curl found. Install one and retry."
    exit 1
fi

echo "Extracting..."
unzip -q diann-linux.zip -d diann-extract

# Find the diann-linux binary (may be in a subdirectory)
DIANN_BIN=$(find diann-extract -name "diann-linux" -type f 2>/dev/null | head -1)
if [ -z "$DIANN_BIN" ]; then
    # Try finding any executable
    DIANN_BIN=$(find diann-extract -type f -executable 2>/dev/null | head -1)
fi

if [ -z "$DIANN_BIN" ]; then
    echo "ERROR: Could not find diann-linux binary in the downloaded archive."
    echo "Contents of archive:"
    find diann-extract -type f
    exit 1
fi

DIANN_DIR=$(dirname "$DIANN_BIN")
echo "Found DIA-NN binary at: $DIANN_BIN"

# Write Dockerfile
cat > Dockerfile << 'DOCKERFILE_EOF'
FROM --platform=linux/amd64 debian:bookworm-slim

# System dependencies for DIA-NN
# libgomp1: OpenMP (parallel processing)
# libstdc++6: C++ standard library
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    libgomp1 \
    libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

# Install .NET SDK 8.0 (required for Thermo .raw file reading)
RUN wget -q https://dot.net/v1/dotnet-install.sh -O /tmp/dotnet-install.sh && \
    chmod +x /tmp/dotnet-install.sh && \
    /tmp/dotnet-install.sh --channel 8.0 --install-dir /usr/share/dotnet && \
    ln -s /usr/share/dotnet/dotnet /usr/bin/dotnet && \
    rm /tmp/dotnet-install.sh

ENV DOTNET_ROOT=/usr/share/dotnet
ENV PATH="$PATH:/usr/share/dotnet"

# Copy DIA-NN binaries
COPY diann-bin/ /opt/diann/

# Make executable and add to PATH
RUN chmod +x /opt/diann/* 2>/dev/null || true && \
    ln -sf /opt/diann/diann-linux /usr/local/bin/diann

# Verify installation
RUN diann --help 2>&1 | head -3 || echo "Note: diann --help returned non-zero (may be normal)"

WORKDIR /work

ENTRYPOINT ["diann"]
DOCKERFILE_EOF

# Copy DIA-NN binaries to build context
mkdir -p diann-bin
cp -r "$DIANN_DIR"/* diann-bin/

echo ""
echo "Building Docker image: ${IMAGE_NAME}"
echo "Platform: linux/amd64 (x86_64)"
echo ""

# Build the image — force amd64 platform
# On Apple Silicon Macs, this builds an x86_64 image that runs via Rosetta 2
docker build --platform linux/amd64 -t "$IMAGE_NAME" .

echo ""
echo "============================================================"
echo "  DIA-NN ${DIANN_VERSION} Docker image built successfully!"
echo "============================================================"
echo ""
echo "  Image name: ${IMAGE_NAME}"
echo ""
echo "  Test with:"
echo "    docker run --rm ${IMAGE_NAME} --help"
echo ""
echo "  To use with DE-LIMP:"
echo "    1. Start DE-LIMP"
echo "    2. Go to 'New Search' tab"
echo "    3. Select 'Local (Docker)' backend"
echo "    4. The image '${IMAGE_NAME}' will be detected automatically"
echo ""
echo "  Optional: Create ~/.delimp_docker.conf to customize:"
echo '    {"diann_image": "'"${IMAGE_NAME}"'", "max_cpus": 32, "max_memory_gb": 128}'
echo ""

# Quick sanity test
echo "Running quick verification..."
if docker run --rm --platform linux/amd64 "$IMAGE_NAME" --help &> /dev/null; then
    echo "  DIA-NN responds to --help"
else
    echo "  Warning: DIA-NN --help returned non-zero exit (may be normal for some versions)"
fi

echo ""
echo "Done!"
