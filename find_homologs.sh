#!/usr/bin/env bash
set -euo pipefail

# Init empty find_homologs.sh
# Init: add script scaffold

[ $# -ge 2 ] || { echo "usage: $0 query.faa subject.fna [outdir]"; exit 1; }
Q=$1; S=$2; O=${3:-homologs_out}; E=${E:-1e-5}; L=${L:-50}

mkdir -p "$O"; DB="$O/db"
makeblastdb -in "$S" -dbtype nucl -out "$DB" >/dev/null 2>&1 || true
tblastn -query "$Q" -db "$DB" -outfmt 6 > "$O/hits.tsv"

awk -v e="$E" -v l="$L" '($9+0)<=e && $4>=l' "$O/hits.tsv" > "$O/hits.filtered.tsv" || true
[ -s "$O/hits.filtered.tsv" ] \
  && cut -f2 "$O/hits.filtered.tsv" | sort -u | wc -l | tr -d ' ' > "$O/hits.count.txt" \
  || echo 0 > "$O/hits.count.txt"

  cat "$O/hits.count.txt" 
