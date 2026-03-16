# Thermo .raw Metadata Parsing — Setup Guide

DE-LIMP can extract instrument metadata from Thermo .raw files via two paths:

## Path 1: Local parsing via `rawrr` (macOS)

For parsing .raw files on your local machine.

### Prerequisites
- Homebrew (`brew`)

### Install

```bash
# 1. Install Mono (.NET runtime for macOS)
brew install mono
```

```r
# 2. Install rawrr from Bioconductor
BiocManager::install("rawrr")

# 3. Download Thermo's RawFileReader assembly (run once)
# This will prompt you to accept Thermo's license agreement
rawrr::installRawFileReaderDLLs()
```

### Verify

```r
rawrr::readFileHeader("/path/to/any/file.raw")
```

---

## Path 2: HPC parsing via ThermoRawFileParser

For parsing remote .raw files over SSH (e.g., on HIVE). DE-LIMP runs this automatically when scanning .raw files via SSH.

### Prerequisites
- Mono on the HPC system (`module load mono` or `apt install mono-complete`)

### Install (pick one)

```bash
# Option A: Conda (recommended)
conda install -c bioconda thermorawfileparser

# Option B: Module (if available on your HPC)
module load thermorawfileparser

# Option C: Manual download
wget https://github.com/compomics/ThermoRawFileParser/releases/latest/download/ThermoRawFileParser.zip
unzip ThermoRawFileParser.zip -d ~/software/thermorawfileparser
chmod +x ~/software/thermorawfileparser/ThermoRawFileParser.sh
echo 'export PATH="$HOME/software/thermorawfileparser:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Verify

```bash
ThermoRawFileParser.sh --version
```

---

## What each tool extracts

| Field | rawrr (local) | ThermoRawFileParser (HPC) |
|-------|:---:|:---:|
| Instrument model | Yes | Yes |
| Serial number | Yes | Yes |
| m/z range | Yes | Yes |
| RT range | Yes | Yes |
| MS1/MS2 counts | Yes | Yes |
| Mass resolution | Yes | Yes |
| Scan filters | Yes | No |
| DIA/DDA detection | Yes (from scan index) | Yes (heuristic) |
| Software version | Yes | Yes |
| LC method/gradient | No | No |

**Note:** Neither tool can extract LC method or gradient information — this is a limitation of Thermo's RawFileReader API. LC conditions must be entered manually for Thermo data.

## How DE-LIMP uses these

- **Local .raw scan**: Uses rawrr if installed, otherwise falls back to extension-based detection ("Orbitrap detected")
- **SSH .raw scan**: Runs ThermoRawFileParser remotely if available, otherwise no metadata (silent fallback)
- **timsTOF .d files**: Work fully without any extra installs (uses built-in SQLite + XML parsing)
