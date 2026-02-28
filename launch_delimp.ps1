# =============================================================================
# DE-LIMP Employee Launcher (Windows PowerShell)
# =============================================================================
# Place this script next to your SSH key and run:
#   .\launch_delimp.ps1
#
# It will:
#   1. Find your SSH key
#   2. SSH to HIVE and install the container if needed
#   3. Clone/update the DE-LIMP repo
#   4. Submit the app via sbatch
#   5. Open an SSH tunnel and launch your browser
# =============================================================================

$ErrorActionPreference = "Stop"

# --- Configuration (edit these if needed) ---
$HIVE_HOST       = "hive.hpc.ucdavis.edu"
$PORT            = 7860
$CORE_DIR        = "/share/genome-center/delimp"
$ACCOUNT         = "genome-center-grp"
$PARTITION       = "high"
$CONFIG_FILE     = ".delimp_config.ps1"
$SETUP_SCRIPT    = "hpc_setup.sh"
$GITHUB_RAW      = "https://raw.githubusercontent.com/bsphinney/DE-LIMP/main"
$MAX_WAIT_NODE   = 600    # seconds to wait for compute node
$MAX_WAIT_APP    = 120    # seconds to wait for app to respond

# --- State ---
$script:TunnelProcess = $null
$script:SlurmJobId    = ""
$script:SshKey        = ""
$script:HiveUser      = ""

# --- Cleanup ---
function Invoke-Cleanup {
    Write-Host ""
    Write-Host "Shutting down..." -ForegroundColor Yellow

    # Cancel SLURM job
    if ($script:SlurmJobId -and $script:SshKey) {
        Write-Host "  Cancelling SLURM job $($script:SlurmJobId)..."
        try {
            & ssh -i $script:SshKey -o StrictHostKeyChecking=accept-new -o ConnectTimeout=5 `
                "$($script:HiveUser)@$HIVE_HOST" `
                "scancel $($script:SlurmJobId) 2>/dev/null; rm -f ~/logs/delimp_node_$($script:SlurmJobId).txt" 2>$null
        } catch {}
    }

    # Kill SSH tunnel
    if ($script:TunnelProcess -and !$script:TunnelProcess.HasExited) {
        Write-Host "  Closing SSH tunnel (PID $($script:TunnelProcess.Id))..."
        try { $script:TunnelProcess.Kill() } catch {}
    }

    Write-Host "Done. Goodbye!" -ForegroundColor Green
}

# Register Ctrl+C handler
[Console]::TreatControlCAsInput = $false
$null = Register-EngineEvent -SourceIdentifier PowerShell.Exiting -Action { Invoke-Cleanup }

function Write-Header {
    Write-Host ""
    Write-Host "============================================" -ForegroundColor Blue
    Write-Host "  DE-LIMP Proteomics - Employee Launcher" -ForegroundColor Blue
    Write-Host "============================================" -ForegroundColor Blue
    Write-Host ""
}

# --- Helper: run command on HIVE ---
function Invoke-HiveSsh {
    param([string]$Command)
    $result = & ssh -i $script:SshKey -o StrictHostKeyChecking=accept-new -o ConnectTimeout=10 `
        "$($script:HiveUser)@$HIVE_HOST" $Command 2>&1
    return ($result -join "`n").Trim()
}

# --- Step 1: Find SSH key ---
function Find-SshKey {
    Write-Host "[1/7] Looking for SSH key..." -ForegroundColor Green

    $scriptDir = Split-Path -Parent $MyInvocation.ScriptName
    if (-not $scriptDir) { $scriptDir = $PWD.Path }

    # Check current directory
    foreach ($name in @("id_ed25519", "id_rsa")) {
        $candidate = Join-Path $scriptDir $name
        if (Test-Path $candidate) {
            $script:SshKey = $candidate
            Write-Host "  Found: $candidate"
            return
        }
    }
    # Check for .pem files in script dir
    $pemFiles = Get-ChildItem -Path $scriptDir -Filter "*.pem" -File -ErrorAction SilentlyContinue
    if ($pemFiles) {
        $script:SshKey = $pemFiles[0].FullName
        Write-Host "  Found: $($script:SshKey)"
        return
    }

    # Fall back to ~/.ssh/
    $sshDir = Join-Path $env:USERPROFILE ".ssh"
    foreach ($name in @("id_ed25519", "id_rsa")) {
        $candidate = Join-Path $sshDir $name
        if (Test-Path $candidate) {
            $script:SshKey = $candidate
            Write-Host "  Found: $candidate"
            return
        }
    }
    $pemFiles = Get-ChildItem -Path $sshDir -Filter "*.pem" -File -ErrorAction SilentlyContinue
    if ($pemFiles) {
        $script:SshKey = $pemFiles[0].FullName
        Write-Host "  Found: $($script:SshKey)"
        return
    }

    Write-Host "No SSH key found!" -ForegroundColor Red
    Write-Host "  Place your SSH key (id_ed25519, id_rsa, or *.pem) next to this script"
    Write-Host "  or in ~/.ssh/"
    exit 1
}

# --- Fix SSH key permissions (Windows requires strict perms) ---
function Repair-SshKeyPermissions {
    $key = $script:SshKey
    if (-not $key -or -not (Test-Path $key)) { return }

    try {
        $acl = Get-Acl $key
        # Check if anyone besides the current user has access
        $otherRules = $acl.Access | Where-Object {
            $_.IdentityReference -notmatch [regex]::Escape($env:USERNAME) -and
            $_.IdentityReference -ne "NT AUTHORITY\SYSTEM" -and
            $_.IdentityReference -ne "BUILTIN\Administrators"
        }
        if ($otherRules) {
            Write-Host "  Fixing SSH key permissions..." -ForegroundColor Yellow
            & icacls $key /inheritance:r /grant:r "${env:USERNAME}:(R)" 2>$null | Out-Null
            Write-Host "  Permissions fixed."
        }
    } catch {
        # Non-fatal — SSH will error if perms are still wrong
        Write-Host "  Warning: Could not check key permissions." -ForegroundColor Yellow
    }
}

# --- Auto-download missing files from GitHub ---
function Get-RequiredFiles {
    $scriptDir = Split-Path -Parent $MyInvocation.ScriptName
    if (-not $scriptDir) { $scriptDir = $PWD.Path }

    $files = @(
        @{ Name = $SETUP_SCRIPT;      Url = "$GITHUB_RAW/$SETUP_SCRIPT" }
        @{ Name = "launch_delimp.ps1"; Url = "$GITHUB_RAW/launch_delimp.ps1" }
    )

    foreach ($f in $files) {
        $local = Join-Path $scriptDir $f.Name
        if (-not (Test-Path $local)) {
            Write-Host "  Downloading $($f.Name) from GitHub..." -ForegroundColor Yellow
            try {
                Invoke-WebRequest -Uri $f.Url -OutFile $local -UseBasicParsing
                Write-Host "  Saved: $local"
            } catch {
                Write-Host "  Failed to download $($f.Name): $($_.Exception.Message)" -ForegroundColor Red
                Write-Host "  Download it manually from: $($f.Url)"
                exit 1
            }
        }
    }
}

# --- Step 2: Get HIVE username ---
function Get-HiveUsername {
    Write-Host "[2/7] Getting HIVE username..." -ForegroundColor Green

    $scriptDir = Split-Path -Parent $MyInvocation.ScriptName
    if (-not $scriptDir) { $scriptDir = $PWD.Path }
    $configPath = Join-Path $scriptDir $CONFIG_FILE

    # Try loading from config
    if (Test-Path $configPath) {
        . $configPath
        if ($SavedHiveUser) {
            $script:HiveUser = $SavedHiveUser
        }
    }

    if ($script:HiveUser) {
        Write-Host "  Using saved username: $($script:HiveUser)"
    } else {
        $script:HiveUser = Read-Host "  Enter your HIVE username"
        if (-not $script:HiveUser) {
            Write-Host "Username cannot be empty." -ForegroundColor Red
            exit 1
        }
        # Save for next time
        "`$SavedHiveUser = '$($script:HiveUser)'" | Out-File -FilePath $configPath -Encoding UTF8
        Write-Host "  Saved to $configPath"
    }
}

# --- Step 3: Check / install container ---
function Test-Container {
    Write-Host "[3/7] Checking container on HIVE..." -ForegroundColor Green

    $hasSif = Invoke-HiveSsh "test -f ~/containers/de-limp.sif && echo yes || echo no"

    if ($hasSif -eq "yes") {
        Write-Host "  Container found."
    } else {
        Write-Host "  Container not found. Installing..." -ForegroundColor Yellow
        Write-Host "  (This downloads ~5 GB and may take 10-20 minutes)"

        # Copy hpc_setup.sh to HIVE
        $scriptDir = Split-Path -Parent $MyInvocation.ScriptName
        if (-not $scriptDir) { $scriptDir = $PWD.Path }
        $setupPath = Join-Path $scriptDir $SETUP_SCRIPT

        if (-not (Test-Path $setupPath)) {
            Write-Host "Cannot find $SETUP_SCRIPT next to this launcher." -ForegroundColor Red
            exit 1
        }

        & scp -i $script:SshKey -o StrictHostKeyChecking=accept-new `
            $setupPath "$($script:HiveUser)@${HIVE_HOST}:~/$SETUP_SCRIPT"

        # Run install (needs TTY for srun)
        & ssh -t -i $script:SshKey -o StrictHostKeyChecking=accept-new `
            "$($script:HiveUser)@$HIVE_HOST" "bash ~/$SETUP_SCRIPT install"

        Write-Host "  Container installed!" -ForegroundColor Green
    }
}

# --- Step 4: Clone/update repo ---
function Update-Repo {
    Write-Host "[4/7] Syncing DE-LIMP repo on HIVE..." -ForegroundColor Green

    $result = Invoke-HiveSsh "if [ -d ~/DE-LIMP/.git ]; then cd ~/DE-LIMP && git pull --ff-only 2>&1 | tail -1; else git clone https://github.com/bsphinney/DE-LIMP.git ~/DE-LIMP 2>&1 | tail -1; fi"
    Write-Host "  $result"
}

# --- Step 5: Check R packages ---
function Test-Packages {
    Write-Host "[5/7] Checking R packages..." -ForegroundColor Green

    $pkgCount = (Invoke-HiveSsh "ls -1d ~/R/delimp-lib/*/ 2>/dev/null | wc -l").Trim()

    if ([int]$pkgCount -lt 3) {
        Write-Host "  Missing R packages. Installing..." -ForegroundColor Yellow

        $scriptDir = Split-Path -Parent $MyInvocation.ScriptName
        if (-not $scriptDir) { $scriptDir = $PWD.Path }
        $setupPath = Join-Path $scriptDir $SETUP_SCRIPT

        & scp -i $script:SshKey -o StrictHostKeyChecking=accept-new `
            $setupPath "$($script:HiveUser)@${HIVE_HOST}:~/$SETUP_SCRIPT" 2>$null

        Invoke-HiveSsh "bash ~/$SETUP_SCRIPT packages"
        Write-Host "  Packages installed!" -ForegroundColor Green
    } else {
        Write-Host "  Found $pkgCount packages - OK."
    }
}

# --- Step 6: Submit via sbatch ---
function Submit-Job {
    Write-Host "[6/7] Submitting DE-LIMP to compute node..." -ForegroundColor Green

    # Ensure hpc_setup.sh is on HIVE
    $scriptDir = Split-Path -Parent $MyInvocation.ScriptName
    if (-not $scriptDir) { $scriptDir = $PWD.Path }
    $setupPath = Join-Path $scriptDir $SETUP_SCRIPT

    & scp -i $script:SshKey -o StrictHostKeyChecking=accept-new `
        $setupPath "$($script:HiveUser)@${HIVE_HOST}:~/$SETUP_SCRIPT" 2>$null

    $submitOutput = & ssh -i $script:SshKey -o StrictHostKeyChecking=accept-new -o ConnectTimeout=10 `
        "$($script:HiveUser)@$HIVE_HOST" "bash ~/$SETUP_SCRIPT sbatch '$CORE_DIR' ~/DE-LIMP 2>&1"
    $submitOutput = ($submitOutput -join "`n").Trim()

    # Parse JOBID:<number>
    if ($submitOutput -match "JOBID:(\d+)") {
        $script:SlurmJobId = $Matches[1]
    } else {
        Write-Host "Failed to submit job. Output:" -ForegroundColor Red
        Write-Host $submitOutput
        Write-Host ""
        Write-Host "Debug: running sbatch command manually..." -ForegroundColor Yellow
        $debugOutput = Invoke-HiveSsh "bash -x ~/$SETUP_SCRIPT sbatch '$CORE_DIR' ~/DE-LIMP 2>&1 | tail -30"
        Write-Host $debugOutput
        exit 1
    }

    Write-Host "  Submitted SLURM job: $($script:SlurmJobId)"

    # Poll for compute node hostname
    Write-Host "  Waiting for compute node allocation..."
    $elapsed = 0
    $node = ""

    while ($elapsed -lt $MAX_WAIT_NODE) {
        $node = (Invoke-HiveSsh "cat ~/logs/delimp_node_$($script:SlurmJobId).txt 2>/dev/null").Trim()
        if ($node) { break }
        Start-Sleep -Seconds 5
        $elapsed += 5
        Write-Host "  Waiting... ($elapsed`s / $MAX_WAIT_NODE`s)" -NoNewline
        Write-Host "`r" -NoNewline
    }
    Write-Host ""

    if (-not $node) {
        Write-Host "Timed out waiting for compute node ($MAX_WAIT_NODE`s)." -ForegroundColor Red
        Write-Host "  Check job status: ssh $($script:HiveUser)@$HIVE_HOST squeue -u $($script:HiveUser)"
        Invoke-Cleanup
        exit 1
    }

    Write-Host "  Running on node: $node" -ForegroundColor Green

    # Open SSH tunnel as background process
    Write-Host "  Opening SSH tunnel (localhost:$PORT -> ${node}:$PORT)..."
    $script:TunnelProcess = Start-Process -PassThru -NoNewWindow ssh `
        -ArgumentList "-i", $script:SshKey, "-o", "StrictHostKeyChecking=accept-new", `
                      "-N", "-L", "${PORT}:${node}:${PORT}", "$($script:HiveUser)@$HIVE_HOST"

    # Poll until app responds
    Write-Host "  Waiting for DE-LIMP to start..."
    $elapsed = 0
    while ($elapsed -lt $MAX_WAIT_APP) {
        try {
            $response = Invoke-WebRequest -Uri "http://localhost:$PORT" -UseBasicParsing -TimeoutSec 3 -ErrorAction SilentlyContinue
            if ($response.StatusCode -eq 200 -or $response.StatusCode -eq 302) { break }
        } catch {}
        Start-Sleep -Seconds 3
        $elapsed += 3
        Write-Host "  Starting... ($elapsed`s / $MAX_WAIT_APP`s)" -NoNewline
        Write-Host "`r" -NoNewline
    }
    Write-Host ""
}

# --- Step 7: Open browser ---
function Open-Browser {
    Write-Host "[7/7] Opening browser..." -ForegroundColor Green

    $url = "http://localhost:$PORT"
    Start-Process $url

    Write-Host ""
    Write-Host "============================================" -ForegroundColor Blue
    Write-Host "  DE-LIMP is running!" -ForegroundColor Blue
    Write-Host "============================================" -ForegroundColor Blue
    Write-Host ""
    Write-Host "  URL:       $url"
    Write-Host "  SLURM Job: $($script:SlurmJobId)"
    Write-Host "  Tunnel:    localhost:$PORT"
    Write-Host ""
    Write-Host "  Press Ctrl+C to stop."
    Write-Host ""

    # Wait forever until Ctrl+C
    try {
        while ($true) {
            Start-Sleep -Seconds 60
        }
    } finally {
        Invoke-Cleanup
    }
}

# --- Main ---
Write-Header
Get-RequiredFiles
Find-SshKey
Repair-SshKeyPermissions
Get-HiveUsername
Test-Container
Update-Repo
Test-Packages
Submit-Job
Open-Browser
