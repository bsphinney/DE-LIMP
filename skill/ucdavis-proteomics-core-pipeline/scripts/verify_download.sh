# =============================================================================
# NO SHEBANG ON PURPOSE: this file is `. sourced` (acquire_tools.sh, setup.sh,
# build_diann_docker.sh, ensure_dotnet8.sh), never executed -- the same reason
# diann_release.sh has none.
#
# verify_download.sh -- ONE definition of "was this download the file the vendor
# published?". Everything the pipeline downloads and then runs goes through it.
#
# The expected sha256 values live in pinned_checksums.txt, beside this file, one line each:
#     <key> <sha256> <where the value came from>
# and each value was copied from the vendor's OWN release metadata (the `digest` field GitHub
# shows for a release asset, or a checksum file the vendor publishes) -- never computed from a
# download we made ourselves, which would only prove the file did not change since we looked.
# A build with no published checksum has NO line there; it is then downloaded and said to be
# "not verified" out loud, never silently trusted.
#
#   sha256_of <file>                      prints the sha256 (sha256sum, shasum or openssl)
#   expected_sha256 <key>                 prints the pinned value, or nothing
#   verify_download <key> <file> <label>  0 = matches the pin
#                                         1 = MISMATCH (a message is printed; the caller deletes the file)
#                                         2 = no pin for this key (a "not verified" warning is printed)
#                                         3 = no tool to compute a sha256 (warning printed)
# =============================================================================

VERIFY_DIR="${VERIFY_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PINNED_CHECKSUMS="${PINNED_CHECKSUMS:-$VERIFY_DIR/pinned_checksums.txt}"

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" 2>/dev/null | cut -d' ' -f1
  elif command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" 2>/dev/null | cut -d' ' -f1
  elif command -v openssl >/dev/null 2>&1; then openssl dgst -sha256 "$1" 2>/dev/null | sed 's/^.*= *//'
  fi
}

expected_sha256() {
  [ -r "$PINNED_CHECKSUMS" ] || return 0
  awk -v k="$1" '$1 == k && $1 !~ /^#/ { print $2; exit }' "$PINNED_CHECKSUMS"
}

verify_download() {
  local key="$1" file="$2" label="${3:-$1}" want got
  want="$(expected_sha256 "$key")"
  if [ -z "$want" ]; then
    echo "  WARNING: $label was NOT verified -- the vendor publishes no checksum we could pin for '$key'. It was downloaded over HTTPS from the vendor's release page and is used as is." >&2
    return 2
  fi
  got="$(sha256_of "$file")"
  if [ -z "$got" ]; then
    echo "  WARNING: $label was NOT verified -- no sha256sum, shasum or openssl on this machine." >&2
    return 3
  fi
  if [ "$got" != "$want" ]; then
    echo "  ERROR: $label does not match the published checksum, so it was NOT used." >&2
    echo "         expected sha256 $want" >&2
    echo "         got      sha256 $got" >&2
    return 1
  fi
  echo "  $label: sha256 verified ($want)" >&2
  return 0
}
