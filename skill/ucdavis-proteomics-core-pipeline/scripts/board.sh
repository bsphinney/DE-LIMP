#!/usr/bin/env bash
# board.sh -- talk to the UC Davis Proteomics Core Project Board (Core staff only).
#
# Runs board_client.py, copied unchanged from the board's own repository (its sha256 is
# pinned in tests/test_board.py), with the Core board's address filled in:
#
#   bash scripts/board.sh connect --owner <their UC Davis email> --wait 600
#   bash scripts/board.sh start-link --prot PROT_0807 --title ... --goal ...
#   bash scripts/board.sh threads
#
# What to do with it, and when: references/board.md. BOARD_URL in the environment wins (a
# test board). BOARD_PYTHON picks the interpreter (default python3; `py -3` on Windows when
# python3 is the Microsoft Store alias).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BOARD_URL="${BOARD_URL:-https://core-board-ucd.azurewebsites.net}"
# shellcheck disable=SC2086  # BOARD_PYTHON may be two words ("py -3")
exec ${BOARD_PYTHON:-python3} "$HERE/board_client.py" "$@"
