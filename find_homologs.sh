#!/usr/bin/env bash
set -euo pipefail

# Init empty find_homologs.sh
# Init: add script scaffold

usage() {
  sed -n '1,/^usage()/{/usage()/!p}; /^# Example:/,$p' "$0" | sed '1d;/^usage()/,$d'
  cat <<'USAGE'
