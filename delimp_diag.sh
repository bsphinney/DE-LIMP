#!/bin/bash
# =============================================================================
# delimp_diag.sh — DE-LIMP Docker install diagnostic collector
# =============================================================================
# Run from inside the DE-LIMP container:
#   docker compose exec delimp bash /tmp/delimp_diag.sh
#
# Or pipe directly from host:
#   docker compose exec delimp bash -c "curl -s <URL_OR_PASTE_SCRIPT> | bash"
#
# Output:
#   /data/diag/delimp_diag_<timestamp>.txt  (visible on Windows host)
#   Also scp'd to hive:/quobyte/proteomics-grp/de-limp/diag/<host>_<ts>.txt
#   if SSH works.
# =============================================================================

set +e  # Don't exit on errors — we want to capture everything

TS=$(date +%Y%m%d_%H%M%S)
HOSTNAME=$(hostname 2>/dev/null || echo unknown)
OUT_DIR="/data/diag"
mkdir -p "$OUT_DIR" 2>/dev/null
OUT="$OUT_DIR/delimp_diag_${HOSTNAME}_${TS}.txt"

# Everything below writes to both stdout (container console) and the log file
exec > >(tee "$OUT") 2>&1

echo "============================================================"
echo "  DE-LIMP Docker Install Diagnostic"
echo "  Timestamp: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "  Container hostname: $HOSTNAME"
echo "  Output: $OUT"
echo "============================================================"

section() { echo ""; echo "### $1 ###"; echo ""; }

section "1. Container identity"
id
echo "cat /etc/os-release:"; cat /etc/os-release 2>/dev/null | head -5
echo "uname -a:"; uname -a
echo ""
echo "Image build metadata (if /.image-build-info exists):"
[ -f /.image-build-info ] && cat /.image-build-info || echo "  (not labeled — image built before build-info metadata was added)"
echo ""
echo "DE-LIMP app version:"
cat /srv/shiny-server/VERSION 2>/dev/null || echo "  (VERSION file not found)"
echo ""
echo "Env vars of interest:"
env | grep -E 'DELIMP|SSH|DOCKER|PATH|HOME|USER' | sort

section "2. Volume mounts visible from inside"
echo "mount | grep -vE 'proc|sysfs|cgroup|tmpfs':"
mount 2>/dev/null | grep -vE 'proc|sysfs|cgroup|tmpfs' | head -30

section "3. Data directory state"
echo "ls -la /data:"
ls -la /data 2>&1 | head -20
echo ""
echo "Sizes:"
du -sh /data/* 2>/dev/null | head -15
echo ""
echo "Disk free:"
df -h /data 2>&1 | head -3
echo ""
echo "Write test on /data/output (used by DIA-NN for predicted library):"
TEST_DIR="/data/output/.delimp_write_test_$$"
if mkdir -p "$TEST_DIR" 2>/dev/null; then
    if echo "test" > "$TEST_DIR/probe.txt" 2>/dev/null && [ -f "$TEST_DIR/probe.txt" ]; then
        echo "  /data/output write: OK"
    else
        echo "  /data/output write: FAILED — DIA-NN speclib save will fail. Check Windows folder permissions + antivirus."
    fi
    rm -rf "$TEST_DIR" 2>/dev/null
else
    echo "  /data/output mkdir: FAILED — parent dir not writable. Check docker-compose volume + Windows permissions."
fi
echo ""
echo "Mount performance hint (9p msize):"
if mount 2>/dev/null | grep -q '/data .*msize='; then
    MSIZE=$(mount 2>/dev/null | grep '/data ' | grep -oE 'msize=[0-9]+' | head -1)
    MSIZE_NUM=$(echo "$MSIZE" | grep -oE '[0-9]+')
    echo "  $MSIZE"
    if [ -n "$MSIZE_NUM" ] && [ "$MSIZE_NUM" -lt 262144 ]; then
        echo "  WARNING: msize < 262144 — 9p bulk I/O will be very slow (~10x slowdown on large spectra files). Consider Docker Desktop WSL2 settings."
    fi
fi

section "4. SSH setup"
echo "data/ssh contents:"
ls -la /data/ssh/ 2>&1 | head -10
echo ""
echo "/tmp/.ssh contents (where entrypoint copies keys + fixes perms):"
ls -la /tmp/.ssh/ 2>&1 | head -10
echo ""
echo "DELIMP_SSH_KEY=$DELIMP_SSH_KEY"
echo "DELIMP_SSH_USER=$DELIMP_SSH_USER"
echo ""
echo "Key file type + line endings (if present):"
for k in /tmp/.ssh/*; do
    [ -f "$k" ] || continue
    echo "  $k:"
    echo "    file type: $(file -b "$k" 2>/dev/null)"
    echo "    perms: $(stat -c '%a %U:%G' "$k" 2>/dev/null)"
    echo "    size: $(stat -c '%s' "$k" 2>/dev/null) bytes"
    # Detect CRLF line endings (Windows corruption that breaks libcrypto)
    if head -c 10000 "$k" 2>/dev/null | grep -q $'\r'; then
        echo "    LINE ENDINGS: CRLF detected (BROKEN). Fix on Windows with:"
        echo "      \$c = [IO.File]::ReadAllText(\"\$PWD\\data\\ssh\\$(basename $k)\") -replace \"\`r\`n\", \"\`n\""
        echo "      [IO.File]::WriteAllText(\"\$PWD\\data\\ssh\\$(basename $k)\", \$c)"
    else
        echo "    LINE ENDINGS: LF (correct)"
    fi
    # First line sanity check
    FIRSTLINE=$(head -1 "$k" 2>/dev/null | tr -d '\r')
    echo "    first line: $FIRSTLINE"
    if [ "$FIRSTLINE" != "-----BEGIN OPENSSH PRIVATE KEY-----" ] && [ "$FIRSTLINE" != "-----BEGIN RSA PRIVATE KEY-----" ] && [ "$FIRSTLINE" != "-----BEGIN EC PRIVATE KEY-----" ]; then
        echo "    WARNING: header doesn't look like a private key. Did you accidentally copy the .pub file?"
    fi
done

section "5. SSH connectivity test to HIVE"
SSH_USER="${DELIMP_SSH_USER:-brettsp}"
SSH_KEY="${DELIMP_SSH_KEY:-/tmp/.ssh/id_ed25519}"
echo "Testing: ssh -i $SSH_KEY ${SSH_USER}@hive.hpc.ucdavis.edu ..."
echo ""
echo "--- verbose output ---"
timeout 15 ssh -v -o BatchMode=yes \
    -o StrictHostKeyChecking=accept-new \
    -o ConnectTimeout=10 \
    -o UserKnownHostsFile=/tmp/.ssh/known_hosts \
    -i "$SSH_KEY" \
    "${SSH_USER}@hive.hpc.ucdavis.edu" \
    "echo SSH_OK; hostname; whoami; groups; pwd; df -h /quobyte/proteomics-grp 2>/dev/null | tail -2" 2>&1 | head -80
SSH_EXIT=$?
echo "--- exit code: $SSH_EXIT ---"

section "6. DIA-NN installation"
echo "which diann: $(which diann 2>&1)"
echo "which diann-linux: $(which diann-linux 2>&1)"
DIANN_BIN=$(which diann 2>/dev/null || which diann-linux 2>/dev/null)
if [ -n "$DIANN_BIN" ]; then
    echo "file: $(file -b "$DIANN_BIN" 2>/dev/null)"
    echo "perms: $(stat -c '%a %U:%G' "$DIANN_BIN" 2>/dev/null)"
    echo "size: $(du -h "$DIANN_BIN" 2>/dev/null | cut -f1)"
    echo ""
    echo "First 3 lines of diann --help:"
    timeout 10 "$DIANN_BIN" --help 2>&1 | head -3
fi
echo ""
echo ".NET runtime:"
which dotnet 2>&1
timeout 5 dotnet --info 2>&1 | head -6

section "7. Library paths (for DIA-NN .raw support)"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "DOTNET_ROOT=$DOTNET_ROOT"
echo ""
echo "Key shared libs in /opt/diann:"
ls -la /opt/diann/ 2>&1 | head -10

section "8. R + processx"
echo "R version:"
R --version 2>&1 | head -2
echo ""
echo "processx version + test:"
R -e '
cat("processx:", as.character(packageVersion("processx")), "\n")
cat("Test processx on /bin/echo...\n")
tryCatch({
  p <- processx::process$new("/bin/echo", c("hello"), stdout = "/tmp/px_test.log")
  p$wait(timeout = 5000)
  cat("  result:", readLines("/tmp/px_test.log"), "\n")
  unlink("/tmp/px_test.log")
}, error = function(e) cat("  FAIL:", e$message, "\n"))
cat("Test processx on DIA-NN...\n")
diann_bin <- Sys.which("diann")
if (nzchar(diann_bin)) {
  tryCatch({
    p <- processx::process$new(diann_bin, "--help", stdout = "/tmp/diann_test.log", stderr = "/tmp/diann_test.log")
    p$wait(timeout = 10000)
    lines <- readLines("/tmp/diann_test.log")
    cat("  first line:", head(lines, 1), "\n")
    cat("  exit:", p$get_exit_status(), "\n")
    unlink("/tmp/diann_test.log")
  }, error = function(e) cat("  FAIL:", e$message, "\n"))
} else cat("  no diann on PATH\n")
' 2>&1 | tail -15

section "9. Network / DNS"
echo "Can resolve hive.hpc.ucdavis.edu:"
getent hosts hive.hpc.ucdavis.edu 2>&1 | head -2
echo ""
echo "Can reach github.com (HTTPS):"
timeout 5 curl -sI https://github.com 2>&1 | head -3
echo ""
echo "Can reach HIVE SSH port (22):"
timeout 5 bash -c 'cat < /dev/tcp/hive.hpc.ucdavis.edu/22' 2>&1 | head -2 | xxd | head -2
echo "  (^ if you see 'SSH-2.0-OpenSSH...' hex above, port 22 is reachable)"

section "10. Recent Shiny / DIA-NN logs"
echo "Last 30 lines of /srv/shiny-server/logs (if any):"
find /srv/shiny-server/logs -type f -name '*.log' -mtime -1 2>/dev/null | head -3 | while read f; do
    echo "--- $f ---"
    tail -30 "$f" 2>/dev/null
done

section "11. Auto-diagnosis — likely issues + fixes"
# Pattern-match the collected log for common failure signatures and emit
# specific fix suggestions. This is the "what do I do now" summary.
echo "Scanning for known failure patterns..."
ISSUES_FOUND=0

# Check 1: CRLF on SSH key
for k in /tmp/.ssh/*; do
    [ -f "$k" ] || continue
    if head -c 10000 "$k" 2>/dev/null | grep -q $'\r'; then
        echo ""
        echo "  [FOUND] SSH key has Windows CRLF line endings — breaks libcrypto."
        echo "     Fix (PowerShell on Windows host):"
        echo "       \$f = \"\$PWD\\data\\ssh\\$(basename $k)\""
        echo "       \$c = [IO.File]::ReadAllText(\$f) -replace \"\`r\`n\", \"\`n\""
        echo "       [IO.File]::WriteAllText(\$f, \$c)"
        echo "       docker compose down; docker compose up -d"
        ISSUES_FOUND=$((ISSUES_FOUND + 1))
    fi
done

# Check 2: SSH permission denied
if grep -q 'Permission denied (publickey' "$OUT" 2>/dev/null; then
    echo ""
    echo "  [FOUND] SSH rejected the key."
    echo "     Check order: (1) public half not on HIVE's ~/.ssh/authorized_keys"
    echo "                 (2) key name in data/ssh/ doesn't match HIVE username"
    echo "                 (3) key passphrase-protected (DE-LIMP can't handle those)"
    echo "     To register the public half on HIVE from Windows PowerShell:"
    echo "       type \$env:USERPROFILE\\.ssh\\<your_key>.pub | ssh brettsp@hive.hpc.ucdavis.edu 'mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys'"
    ISSUES_FOUND=$((ISSUES_FOUND + 1))
fi

# Check 3: DIA-NN binary not found
if ! [ -x "$(which diann 2>/dev/null)" ] && ! [ -x "$(which diann-linux 2>/dev/null)" ]; then
    echo ""
    echo "  [FOUND] DIA-NN binary missing from container."
    echo "     Fix: rebuild on Windows:"
    echo "       .\\build_diann_docker.ps1"
    echo "       docker compose down; docker compose up -d --build"
    ISSUES_FOUND=$((ISSUES_FOUND + 1))
fi

# Check 4: /data/output not writable
TEST_W="/data/output/.w_$$"
if ! mkdir -p "/data/output" 2>/dev/null || ! echo t > "$TEST_W" 2>/dev/null; then
    echo ""
    echo "  [FOUND] /data/output not writable — DIA-NN outputs and predicted library saves will fail."
    echo "     Check: Windows folder permissions on <DE-LIMP>/data/, and whether antivirus is blocking container writes."
    echo "     Verify from PowerShell: echo test > \"\$PWD\\data\\output\\test.txt\""
    ISSUES_FOUND=$((ISSUES_FOUND + 1))
fi
rm -f "$TEST_W" 2>/dev/null

# Check 5: Port 22 unreachable
if ! timeout 3 bash -c '</dev/tcp/hive.hpc.ucdavis.edu/22' 2>/dev/null; then
    echo ""
    echo "  [FOUND] Cannot reach hive.hpc.ucdavis.edu on port 22."
    echo "     Your network (corporate firewall, VPN) is blocking outbound SSH."
    echo "     Options:"
    echo "       (1) Connect to UC Davis VPN"
    echo "       (2) Ask IT to allow outbound 22 to hive.hpc.ucdavis.edu (169.237.253.41)"
    echo "       (3) Use Local DIA-NN backend instead — doesn't need HPC connection"
    ISSUES_FOUND=$((ISSUES_FOUND + 1))
fi

# Check 6: Low 9p msize (performance, not correctness)
if mount 2>/dev/null | grep -q '/data .*msize='; then
    MSIZE_NUM=$(mount 2>/dev/null | grep '/data ' | grep -oE 'msize=[0-9]+' | head -1 | grep -oE '[0-9]+')
    if [ -n "$MSIZE_NUM" ] && [ "$MSIZE_NUM" -lt 262144 ]; then
        echo ""
        echo "  [INFO] 9p msize=$MSIZE_NUM is small — large file I/O will be slow."
        echo "     Not a correctness issue, but DIA-NN speclib writes and raw file reads will be ~10x slower than they could be."
        echo "     Fix: set /etc/wsl.conf in your WSL2 distro (needs WSL restart):"
        echo "       [automount]"
        echo "       options = \"metadata,msize=262144,cache=mmap\""
        ISSUES_FOUND=$((ISSUES_FOUND + 1))
    fi
fi

# Check 7: SSH key missing / empty
if [ -z "$(ls /tmp/.ssh/ 2>/dev/null | grep -v gitkeep)" ]; then
    echo ""
    echo "  [FOUND] No SSH key in container. HPC submission will not work."
    echo "     Fix:"
    echo "       1. Place your HIVE private key at <DE-LIMP>\\data\\ssh\\<hive_username>"
    echo "       2. docker compose down; docker compose up -d"
    ISSUES_FOUND=$((ISSUES_FOUND + 1))
fi

if [ "$ISSUES_FOUND" -eq 0 ]; then
    echo ""
    echo "  No known failure patterns detected. If the app still misbehaves, share this full log."
else
    echo ""
    echo "  → $ISSUES_FOUND issue(s) flagged above. Work through the fixes in order, then rerun this diagnostic."
fi

echo ""
echo "============================================================"
echo "  Diagnostic complete"
echo "  Saved to: $OUT"
echo "============================================================"

# Attempt upload to HIVE if SSH worked (exit 0 from section 5)
if [ "$SSH_EXIT" -eq 0 ] && [ -n "$SSH_KEY" ] && [ -f "$SSH_KEY" ]; then
    REMOTE_DIR="/quobyte/proteomics-grp/de-limp/diag"
    echo ""
    echo "Uploading diagnostic to HIVE: ${SSH_USER}@hive:${REMOTE_DIR}/"
    timeout 30 ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new \
        -o UserKnownHostsFile=/tmp/.ssh/known_hosts \
        -i "$SSH_KEY" "${SSH_USER}@hive.hpc.ucdavis.edu" \
        "mkdir -p $REMOTE_DIR" 2>&1 | head -3
    timeout 60 scp -o BatchMode=yes -o StrictHostKeyChecking=accept-new \
        -o UserKnownHostsFile=/tmp/.ssh/known_hosts \
        -i "$SSH_KEY" "$OUT" \
        "${SSH_USER}@hive.hpc.ucdavis.edu:${REMOTE_DIR}/" 2>&1 | head -5
    SCP_EXIT=$?
    if [ "$SCP_EXIT" -eq 0 ]; then
        echo "Upload succeeded: ${REMOTE_DIR}/$(basename $OUT)"
    else
        echo "Upload FAILED (exit $SCP_EXIT). The log is still available at $OUT on the Windows host."
    fi
else
    echo ""
    echo "SSH failed in section 5 — not uploading. Please share the log manually."
    echo "On the Windows host, the file is at: <DE-LIMP-dir>\\data\\diag\\$(basename $OUT)"
fi
