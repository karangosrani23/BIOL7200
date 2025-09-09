#!/usr/bin/env bash
set -euo pipefail
usage(){ echo "Usage: $0 -q QUERY.faa -s SUBJECT.fna -o out.txt [-e 1e-3] [-L 30] [-t 2]"; exit 1; }

E=1e-3
L=30
T=2
Q="" S="" OUT=""

while getopts ":q:s:o:e:L:t:h" opt; do
  case $opt in
    q) Q=$OPTARG ;;
    s) S=$OPTARG ;;
    o) OUT=$OPTARG ;;
    e) E=$OPTARG ;;
    L) L=$OPTARG ;;
    t) T=$OPTARG ;;
    h) usage ;;
    \?|:) usage ;;
  esac
done
shift $((OPTIND-1))
if [ -z "$Q" ] && [ $# -ge 3 ]; then Q=$1; S=$2; OUT=$3; fi
[ -n "$Q" ] && [ -n "$S" ] && [ -n "$OUT" ] || usage
[ -f "$Q" ] || { echo "Query not found: $Q" >&2; exit 1; }
[ -f "$S" ] || { echo "Subject not found: $S" >&2; exit 1; }
command -v makeblastdb >/dev/null 2>&1 || { echo "makeblastdb not found" >&2; exit 1; }
command -v tblastn     >/dev/null 2>&1 || { echo "tblastn not found"     >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
db="$tmp/db"
makeblastdb -in "$S" -dbtype nucl -out "$db" >/dev/null
tblastn -query "$Q" -db "$db" -num_threads "$T" -outfmt '6 length evalue sseqid' > "$tmp/hits.tsv" || true

awk -v e="$E" -v l="$L" '($2+0)<=e && $1>=l{print $3}' "$tmp/hits.tsv" | sort -u > "$OUT"
count=$(wc -l < "$OUT" | tr -d " ")
echo "$count"
