# Reading Thermo `.raw` on WSL2 / Linux — Lessons Learned

**Audience:** Claude/sibling tools (STAN, etc.) needing to read Thermo `.raw` mass-spec files on Linux or WSL2 Ubuntu.
**Source:** DE-LIMP v3.10.18–30 install-stack triage on 2026-05-07. 27 hotfix versions to land at the actual root cause; this doc lets the next tool skip 26 of them.
**Maintained by:** Brett Phinney's DE-LIMP project — see `delimp_wsl_setup.sh` for the canonical implementation.

---

## TL;DR — what to install

For any tool whose Thermo `.raw` reader is `RawFileReader.dll` (DIA-NN, ProteoWizard's msconvert in some configs, or anything wrapping Thermo's RawFileReader assembly):

1. **Install the .NET 8 *SDK*** (NOT just the runtime — see §2)
2. **Install .NET 8 system dependencies** explicitly via apt (NOT pulled by `dotnet-install.sh` — see §3)
3. **Set `LD_LIBRARY_PATH`** to your tool's directory at run time (the C++ binary loads bundled libs from there)
4. **Set `DOTNET_ROOT=/usr/share/dotnet`** so the .NET runtime resolver finds the SDK
5. **Verify with `ldd | grep "not found"`** — empty output = libraries resolve

If you do all five, Thermo `.raw` reading on Linux works.

---

## 1. The error you're chasing

Symptom: `cannot read .raw files, please download and install .NET Runtime .NET SDK 8.0.407 or later`

Or the cascade of:
- `Loading run X.raw`
- `No MS2 spectra: aborting`
- `ERROR: cannot load the file, skipping`
- `WARNING: not enough peptides for normalisation`
- Search produces 0 precursors / empty library

These are downstream of the same root cause: the .raw reader can't load the Thermo binary format. Almost always a .NET install problem.

## 2. The .NET 8 SDK requirement (NOT runtime)

DIA-NN 2.x explicitly says **"SDK 8.0.407 or later"**. Even though SDK and runtime sound interchangeable, the tool's check looks for SDK presence. Common mistakes:

- ❌ `apt-get install dotnet-runtime-8.0` — runtime only, not SDK
- ❌ `dotnet-install.sh --runtime dotnet --channel 8.0` — runtime only
- ❌ Using `mcr.microsoft.com/dotnet/runtime:8.0` as a Docker base image
- ✅ `apt-get install dotnet-sdk-8.0` (when available — see §4)
- ✅ `dotnet-install.sh --channel 8.0` (no `--runtime` flag → installs SDK)
- ✅ `mcr.microsoft.com/dotnet/sdk:8.0` Docker base image (~700 MB; runtime is ~200 MB)

Verification:
```bash
dotnet --list-runtimes  # both runtime + SDK installs show this
dotnet --list-sdks      # ONLY shows entries if SDK installed
```

If `--list-runtimes` shows `Microsoft.NETCore.App 8.0.x` but `--list-sdks` is empty → you have runtime, not SDK. The Thermo reader will fail.

## 3. .NET 8 system dependencies (apt deps NOT pulled by dotnet-install.sh)

`dotnet-install.sh` explicitly warns *"the script does not resolve dependencies during installation"*. Without these, `dotnet --list-runtimes` works (just file headers) but **executing a .NET binary fails silently with no output**:

```bash
sudo apt-get install -y --no-install-recommends \
    libc6 libgcc-s1 libssl3 libstdc++6 libunwind8 zlib1g \
    libgssapi-krb5-2 liblttng-ust1
```

Plus libicu, whose package name varies by Ubuntu release:
- Ubuntu 26.04 → `libicu76`
- Ubuntu 24.04 → `libicu74`
- Ubuntu 22.04 → `libicu70`

Probe in order:
```bash
for ic in libicu76 libicu74 libicu72 libicu71 libicu70; do
    if apt-cache show "${ic}" >/dev/null 2>&1; then
        sudo apt-get install -y --no-install-recommends "${ic}"
        break
    fi
done
```

## 4. The Microsoft apt-repo trap on new Ubuntu releases

When Ubuntu ships a new release (e.g. 26.04 in April 2026), Microsoft's apt channel for that release **doesn't immediately have .NET 8 packages**. The script can:

1. Successfully add `https://packages.microsoft.com/config/ubuntu/26.04/packages-microsoft-prod.deb` (the GPG/repo file)
2. `apt-get update` succeeds
3. `apt-get install dotnet-sdk-8.0` fails: `Unable to locate package`

This is the killer for fresh Windows 11/2026 boxes running fresh WSL Ubuntu. **Always have a fallback to `dotnet-install.sh`.**

## 5. The 4-tier install (canonical from DE-LIMP `delimp_wsl_setup.sh`)

```bash
install_dotnet8_sdk() {
    # Tier 1: detect existing 8.x SDK
    if command -v dotnet >/dev/null 2>&1; then
        if dotnet --list-sdks 2>/dev/null | grep -qE '^8\.'; then
            return 0
        fi
        # runtime present but no SDK -> proceed
    fi
    
    # Tier 2: try apt with multiple package names (naming has shifted)
    for pkg in dotnet-sdk-8.0 dotnet-sdk-8 dotnet8-sdk; do
        if apt-cache show "${pkg}" >/dev/null 2>&1; then
            if sudo apt-get install -y "${pkg}"; then return 0; fi
        fi
    done
    
    # Tier 3: Microsoft's apt repo + sweep again
    sudo apt-get install -y wget
    local urel="$(lsb_release -rs)"
    wget -q "https://packages.microsoft.com/config/ubuntu/${urel}/packages-microsoft-prod.deb" \
        -O /tmp/packages-microsoft-prod.deb || true
    sudo dpkg -i /tmp/packages-microsoft-prod.deb >/dev/null 2>&1 || true
    sudo apt-get update -qq || true
    for pkg in dotnet-sdk-8.0 dotnet-sdk-8 dotnet8-sdk; do
        if apt-cache show "${pkg}" >/dev/null 2>&1; then
            if sudo apt-get install -y "${pkg}"; then return 0; fi
        fi
    done
    
    # Tier 4: Microsoft's official dotnet-install.sh
    # NOTE: --channel 8.0 (NOT --version 8.0; that's interpreted as
    # the literal version "8.0" and 404s)
    # NOTE: NO --runtime flag (that installs runtime only; we want SDK)
    curl -sSL https://dot.net/v1/dotnet-install.sh -o /tmp/dotnet-install.sh
    chmod +x /tmp/dotnet-install.sh
    sudo /tmp/dotnet-install.sh --channel 8.0 --install-dir /usr/share/dotnet
    sudo ln -sf /usr/share/dotnet/dotnet /usr/local/bin/dotnet
    
    # CRITICAL: install system deps now (dotnet-install.sh doesn't)
    sudo apt-get install -y --no-install-recommends \
        libc6 libgcc-s1 libssl3 libstdc++6 libunwind8 zlib1g \
        libgssapi-krb5-2 liblttng-ust1
    for ic in libicu76 libicu74 libicu72 libicu71 libicu70; do
        if apt-cache show "${ic}" >/dev/null 2>&1; then
            sudo apt-get install -y --no-install-recommends "${ic}"
            break
        fi
    done
    
    # Verify
    if dotnet --list-sdks 2>/dev/null | grep -qE '^8\.'; then
        return 0
    fi
    return 1
}
```

## 6. Verification (always run after install)

Five checks. If any fail, the tool will fail at runtime with cryptic errors.

```bash
verify_thermo_raw_runtime() {
    local TOOL_DIR="$1"   # directory with the .raw-reading tool's binary
    local TOOL_BIN="$2"   # binary name relative to TOOL_DIR
    
    # 1. dotnet 8 SDK present
    dotnet --list-sdks 2>/dev/null | grep -qE '^8\.' \
        || { echo "FAIL: no .NET 8 SDK"; return 1; }
    
    # 2. tool binary exists and is executable
    [ -x "${TOOL_DIR}/${TOOL_BIN}" ] \
        || { echo "FAIL: ${TOOL_BIN} missing or not executable"; return 1; }
    
    # 3. RawFileReader DLLs bundled
    n=$(find "${TOOL_DIR}" -maxdepth 2 -name '*RawFileReader*' 2>/dev/null | wc -l)
    [ "${n}" -ge 1 ] \
        || { echo "FAIL: no RawFileReader DLLs in ${TOOL_DIR}"; return 1; }
    
    # 4. all dynamic libs resolve
    missing=$(ldd "${TOOL_DIR}/${TOOL_BIN}" 2>&1 | grep 'not found' || true)
    [ -z "${missing}" ] \
        || { echo "FAIL: ldd reports missing libs:"; echo "${missing}"; return 1; }
    
    # 5. binary executes successfully (with proper env)
    LD_LIBRARY_PATH="${TOOL_DIR}:${LD_LIBRARY_PATH}" \
        DOTNET_ROOT=/usr/share/dotnet \
        "${TOOL_DIR}/${TOOL_BIN}" --help >/dev/null 2>&1 \
        || { echo "FAIL: ${TOOL_BIN} --help exited non-zero"; return 1; }
    
    echo "OK: Thermo .raw runtime verified"
}
```

## 7. WSL2-specific gotchas

### 7.1 9P filesystem performance for binary reads
Files at `/mnt/c/...` (Windows side) are accessed through the 9P protocol. **Random-access reads are slow and occasionally unreliable** for the seek-heavy I/O Thermo's RawFileReader does. Symptoms:
- Search succeeds on small files but hangs/fails on large ones
- Inconsistent failures across runs of the same data

**Workaround:** copy or stage `.raw` files onto the WSL native ext4 filesystem (e.g. `~/raw_data/`) before running. Disk space cost trades against reliability/speed.

### 7.2 Permission/ownership
Files copied from Windows to WSL via `/mnt/c/` come in with `root:root` ownership and weird permission bits. If your tool runs as a non-root user, may fail to open. Fix: `chown -R $USER:$USER ~/raw_data/` after copy.

### 7.3 Modern WSL bootstrap
On a fresh Windows 11 box, WSL2 is NOT installed by default. Need:
```powershell
# Run as Administrator
wsl.exe --install --no-distribution     # adds Virtual Machine Platform feature + WSL2 kernel
# REBOOT
wsl --install -d Ubuntu                   # creates the Ubuntu distro
```
And BIOS-side virtualization (Intel VT-x / AMD-V / "SVM Mode") must be enabled. Check via `systeminfo | findstr "Virtualization"` — must say `Yes`.

### 7.4 Probing whether Ubuntu actually exists
The naive check `wsl -d Ubuntu -e true` returns negative-ish exit codes for `WSL_E_DISTRO_NOT_FOUND` that Windows batch's `if errorlevel 1` interprets as success. Use a sentinel-string probe instead:
```batch
wsl -d Ubuntu -e bash -c "echo __SENTINEL_OK__" 2>&1 | findstr /c:"__SENTINEL_OK__" >nul
if errorlevel 1 (
    echo Ubuntu not installed
)
```

## 8. ThermoRawFileParser specifically

If you're using ProteoWizard's [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) (separate from DIA-NN's bundled reader):

- Available as a [BioContainer](https://quay.io/repository/biocontainers/thermorawfileparser) — the easiest path: `docker run --rm -v $PWD:/data quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0 ThermoRawFileParser.sh -i /data/test.raw`
- Native install: same .NET 8 SDK requirement; same ldd checks apply
- Output formats: `.mzML`, `.mgf`, `.parquet` (newer versions) — pick by `-f` flag

For DE-LIMP's TIC extraction over SSH, we use ThermoRawFileParser's BioContainer with `-m=0 -f=0` to extract metadata as JSON without converting the full file.

## 9. Common failure modes — quick lookup

| Symptom | Cause | Fix |
|---|---|---|
| "cannot read .raw files, please install .NET SDK 8.0.407+" | Runtime installed, not SDK | §2 — install SDK |
| `dotnet --list-sdks` empty but `--list-runtimes` populated | Runtime-only install | §4 tier 4 without `--runtime` flag |
| Binary exits silently (no stdout, no stderr, exit 0 or 134) | Missing system libs | §3 — apt-install libicu / libssl3 / etc. |
| `ldd | grep 'not found'` shows libicu, libssl, ... | Same as above | §3 |
| "No MS2 spectra: aborting" on every .raw | Reader can't decode the binary | Re-verify SDK + system deps (§6) |
| Search hangs on large `.raw` but works on small ones | WSL2 9P filesystem | §7.1 — copy to native ext4 |
| `dotnet-install.sh: Could not find with version = 8.0` | Wrong flag | Use `--channel 8.0` not `--version 8.0` |
| `wsl --install -d Ubuntu`: `HCS_E_HYPERV_NOT_INSTALLED` | BIOS virtualization off OR Virtual Machine Platform feature missing | `wsl --install --no-distribution` (admin) + reboot; verify BIOS VT-x |
| Launcher silent-no-ops | `wsl -d Ubuntu -e true` returns weirdly | §7.4 — sentinel-string probe |

## 10. Cross-references in DE-LIMP source

For the canonical implementation, see:
- `delimp_wsl_setup.sh:install_dotnet8_runtime()` — 4-tier install
- `delimp_wsl_setup.sh:install_dotnet_system_deps()` — apt deps
- `delimp_wsl_setup.sh:verify_diann_runtime()` — 5-check verification
- `Launch_DE-LIMP_WSL.bat` — sentinel probe for Ubuntu existence
- `build_diann_docker.sh` — switched from `runtime:8.0` to `sdk:8.0` base image
- `WINDOWS_WSL_INSTALL.md` — user-facing troubleshooting guide

## 11. Decision log

- **2026-05-07** — DE-LIMP shipped 27 hotfix versions (3.10.4 → 3.10.30) tracking down what turned out to be a single root cause: DIA-NN 2.x requires the .NET 8 SDK, not just the runtime. Each layer of fix surfaced the next problem because the earlier layers had been silently masking failures. This doc captures everything learned so the next tool can skip the layers.
