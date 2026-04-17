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
echo "Key file type (if present):"
for k in /tmp/.ssh/*; do
    [ -f "$k" ] || continue
    echo "  $k: $(file -b "$k" 2>/dev/null)"
    echo "  perms: $(stat -c '%a %U:%G' "$k" 2>/dev/null)"
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
