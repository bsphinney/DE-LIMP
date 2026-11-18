# =============================================================================
# NO SHEBANG ON PURPOSE: sourced by detect_env.sh and check_access.sh, never executed.
#
# hive_marker.sh -- ONE shell definition of "is this machine HIVE?" and "is SLURM here?".
# (run_search.cluster_is_hive() is the same rule for Python; keep the two in step.)
#
# Mode used to be decided by `sbatch` being on PATH. That is wrong both ways:
#   * a non-login shell on HIVE (an agent's `ssh hive 'cmd'`, a shell after `module purge`) has NO
#     sbatch on PATH, so a HIVE node looked like "no cluster";
#   * any other SLURM cluster, or a Mac with the SLURM client, looked like HIVE, and was then given
#     HIVE's account and queue names.
# HPC@UCD names its cluster in /opt/hpccf/etc/cluster_name.conf ("hive", "farm", "franklin"). That
# file is HIVE's own marker. Where it is missing (a machine HPCCF does not manage) the Core's group
# folder is the fallback marker, as detect_env.sh has always used.
#
#   slurm_present   sbatch on PATH, or at HPC@UCD's fixed path (HPCCF_SLURM_BIN overrides it)
#   cluster_name    the name in the cluster file, lower case, or nothing
#   is_hive         0 when this is HIVE
# Test overrides: UCD_HPC_CLUSTER_CONF, HPCCF_SLURM_BIN, HIVE_GROUP_DIR.
# =============================================================================

_hm_conf()     { printf '%s' "${UCD_HPC_CLUSTER_CONF:-/opt/hpccf/etc/cluster_name.conf}"; }
_hm_slurmbin() { printf '%s' "${HPCCF_SLURM_BIN:-/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin}"; }

slurm_present() {
  command -v sbatch >/dev/null 2>&1 || [ -x "$(_hm_slurmbin)/sbatch" ]
}

# The name, by ONE rule (run_search.cluster_name_from_conf is the same; tests/
# test_hive_marker_agree_2151.py runs both on the same files): HPCCF's own ucd-hpc skill reads the
# file as shell (`. cluster_name.conf; echo $CLUSTER_NAME`), so a `CLUSTER_NAME=<name>` line
# (optionally `export`ed, optionally quoted) gives the name; a file with no such line is the bare
# name -- its first line that is not blank and not a # comment, with its white space removed.
# Lower case either way. (2.15.1 review: the shell read 40 characters of the WHOLE file and Python
# 200, so a longer header could make them disagree, and a comment naming "hive" made any cluster
# HIVE.)
cluster_name() {
  # sed is the ONE outside program here (no grep, tr or cut), and the function never fails:
  # detect_env.sh runs under `set -e`, and the 2.15.1 HIVE suite found `cut` missing from a stub
  # PATH ended the script before it printed any JSON (exit 127, empty stdout).
  local f n; f="$(_hm_conf)"
  [ -r "$f" ] || return 0
  # "=<name>" from the first CLUSTER_NAME line ("=" when its value is empty); nothing when no line
  n="$(sed -n "/^[[:space:]]*\(export[[:space:]]\{1,\}\)\{0,1\}CLUSTER_NAME[[:space:]]*=/{s/^[[:space:]]*\(export[[:space:]]\{1,\}\)\{0,1\}CLUSTER_NAME[[:space:]]*=[[:space:]]*[\"']\{0,1\}\([^\"'[:space:]#]*\).*/=\2/;p;q;}" "$f" 2>/dev/null)" || n=""
  if [ -n "$n" ]; then
    n="${n#=}"
  else
    n="$(sed -n '/^[[:space:]]*#/d;/^[[:space:]]*$/d;s/[[:space:]]//g;p;q' "$f" 2>/dev/null)" || n=""
  fi
  printf '%s\n' "$n" | sed 'y/ABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyz/' 2>/dev/null || true
}

is_hive() {
  local n; n="$(cluster_name)"
  if [ -n "$n" ]; then case "$n" in *hive*) return 0 ;; *) return 1 ;; esac; fi
  # no cluster file: the Core's group folder, readable
  [ -d "${HIVE_GROUP_DIR:-/quobyte/proteomics-grp}" ]
}
