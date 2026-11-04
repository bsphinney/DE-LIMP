#!/usr/bin/env bash
# board.sh -- talk to the UC Davis Proteomics Core Project Board (Core staff only).
#
# Runs board_client.py, copied unchanged from the board's own repository (its sha256 is
# pinned in tests/test_board.py), with the Core board's address filled in:
#
#   bash scripts/board.sh connect --owner <their HIVE login>@ucdavis.edu --wait 600
#   bash scripts/board.sh start-link --prot PROT_0807 --title ... --goal ...
#   bash scripts/board.sh threads
#
# What to do with it, and when: references/board.md. BOARD_URL in the environment wins (a
# test board). BOARD_PYTHON picks the interpreter (default python3; `py -3` on Windows when
# python3 is the Microsoft Store alias).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BOARD_URL="${BOARD_URL:-https://core-board-ucd.azurewebsites.net}"
# 2.14.2: the "delivered" and "finished" rows (references/board.md section 3) say the results are
# ready -- `where ... --session-folder S --report ...` (or --bioshare-url). A WARNING, never a
# refusal, when that session, here, is not finalized (no output/FINALIZED.json).
if [ "${1:-}" = "where" ]; then
  sf="" ready="" prev=""
  for arg in "$@"; do
    [ "$prev" = "--session-folder" ] && sf="$arg"
    case "$arg" in
      --report|--report=*|--bioshare-url|--bioshare-url=*) ready=1 ;;
      --session-folder=*) sf="${arg#--session-folder=}" ;;
    esac
    prev="$arg"
  done
  if [ -n "$sf" ] && [ -n "$ready" ] && [ -d "$sf/output" ] && [ ! -f "$sf/output/FINALIZED.json" ]; then
    echo "[board] WARNING: this session is not finalized (no output/FINALIZED.json in $sf): run" \
         "session.py finalize before saying it is delivered or finished (AGENTS.md, README," \
         "Methods and the bundle come from it)" >&2
  fi
fi
# shellcheck disable=SC2086  # BOARD_PYTHON may be two words ("py -3")
exec ${BOARD_PYTHON:-python3} "$HERE/board_client.py" "$@"
