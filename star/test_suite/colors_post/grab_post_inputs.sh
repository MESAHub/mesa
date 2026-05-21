#!/usr/bin/env bash
set -euo pipefail

src="${1:-}"

if [[ -z "$src" ]]; then
    echo "usage: $0 /path/to/source/test_case"
    exit 2
fi

if [[ ! -d "$src" ]]; then
    echo "error: source directory does not exist: $src"
    exit 1
fi

src="$(realpath "$src")"
dst="$(pwd)"
stamp="$(date +%Y%m%d_%H%M%S)"

echo "source: $src"
echo "dest  : $dst"
echo

# ------------------------------------------------------------
# Find inlist
# Prefer literal ./inlist, otherwise accept exactly one inlist*
# ------------------------------------------------------------
if [[ -f "$src/inlist" ]]; then
    inlist_src="$src/inlist"
else
    mapfile -t inlist_candidates < <(
        find "$src" -maxdepth 1 -type f -name 'inlist*' | sort
    )

    if (( ${#inlist_candidates[@]} == 0 )); then
        echo "error: no inlist or inlist* found in $src"
        exit 1
    elif (( ${#inlist_candidates[@]} > 1 )); then
        echo "error: multiple inlist candidates found:"
        printf '  %s\n' "${inlist_candidates[@]}"
        echo
        echo "Refusing to guess. Keep one, or copy manually."
        exit 1
    else
        inlist_src="${inlist_candidates[0]}"
    fi
fi

# ------------------------------------------------------------
# Find history.data
# Prefer LOGS/history.data, otherwise accept exactly one match
# ------------------------------------------------------------
if [[ -f "$src/LOGS/history.data" ]]; then
    hist_src="$src/LOGS/history.data"
else
    mapfile -t hist_candidates < <(
        find "$src" -maxdepth 4 -type f -name 'history.data' | sort
    )

    if (( ${#hist_candidates[@]} == 0 )); then
        echo "error: no history.data found under $src"
        exit 1
    elif (( ${#hist_candidates[@]} > 1 )); then
        echo "error: multiple history.data candidates found:"
        printf '  %s\n' "${hist_candidates[@]}"
        echo
        echo "Refusing to guess. Copy manually or point at a more specific folder."
        exit 1
    else
        hist_src="${hist_candidates[0]}"
    fi
fi

echo "using inlist : $inlist_src"
echo "using history: $hist_src"
echo

# ------------------------------------------------------------
# Check required history columns
# MESA history column names are normally on line 6
# ------------------------------------------------------------
header="$(awk 'NR==6 {print; exit}' "$hist_src")"

missing=0

if ! grep -qw 'surface_h1' <<< "$header"; then
    echo "error: history file is missing column: surface_h1"
    missing=1
fi

if ! grep -qw 'surface_he4' <<< "$header"; then
    echo "error: history file is missing column: surface_he4"
    missing=1
fi

if (( missing != 0 )); then
    echo
    echo "This history.data was not produced with the required columns."
    echo "Add surface_h1 and surface_he4 to the producer's history_columns.list,"
    echo "rerun that model, then run this script again."
    exit 1
fi

# ------------------------------------------------------------
# Backup existing destination files
# ------------------------------------------------------------
mkdir -p "$dst/LOGS"

if [[ -f "$dst/inlist" ]]; then
    cp -p "$dst/inlist" "$dst/inlist.backup_$stamp"
    echo "backup: ./inlist -> ./inlist.backup_$stamp"
fi

if [[ -f "$dst/LOGS/history.data" ]]; then
    cp -p "$dst/LOGS/history.data" "$dst/LOGS/history.data.backup_$stamp"
    echo "backup: ./LOGS/history.data -> ./LOGS/history.data.backup_$stamp"
fi

# ------------------------------------------------------------
# Copy into canonical post-processing locations
# ------------------------------------------------------------
cp -p "$inlist_src" "$dst/inlist"
cp -p "$hist_src" "$dst/LOGS/history.data"

echo
echo "copied:"
echo "  $inlist_src"
echo "    -> ./inlist"
echo "  $hist_src"
echo "    -> ./LOGS/history.data"
echo
echo "done"
