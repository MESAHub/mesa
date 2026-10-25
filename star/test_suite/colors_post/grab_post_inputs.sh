#!/usr/bin/env bash
set -euo pipefail

# copy the inlists and history file of a finished MESA run into this
# directory, so that ./rn post-processes that run.
#
# all inlist* files are copied, since a run's inlist usually reads its
# &colors namelist from another inlist file by a relative name.

src="${1:-}"

if [[ -z "$src" ]]; then
    echo "usage: $0 /path/to/source/run"
    exit 2
fi

if [[ ! -d "$src" ]]; then
    echo "error: source directory does not exist: $src"
    exit 1
fi

src="$(cd "$src" && pwd)"
dst="$(pwd)"
stamp="$(date +%Y%m%d_%H%M%S)"

echo "source: $src"
echo "dest  : $dst"
echo

if [[ ! -f "$src/inlist" ]]; then
    echo "error: no inlist found in $src"
    exit 1
fi

hist_src="$src/LOGS/history.data"
if [[ ! -f "$hist_src" ]]; then
    echo "error: no LOGS/history.data found in $src"
    echo "Copy the history file manually, or set history_file in &post_proc."
    exit 1
fi

# ------------------------------------------------------------
# Check required history columns
# MESA history column names are on line 6
# ------------------------------------------------------------
header="$(awk 'NR==6 {print; exit}' "$hist_src")"

missing=0
for col in surface_h1 surface_he4; do
    if ! grep -qw "$col" <<< "$header"; then
        echo "error: history file is missing column: $col"
        missing=1
    fi
done

if (( missing != 0 )); then
    echo
    echo "This history.data was not produced with the required columns."
    echo "Add surface h1 and surface he4 to the run's history_columns.list,"
    echo "rerun that model, then run this script again."
    exit 1
fi

# ------------------------------------------------------------
# Back up and copy
# ------------------------------------------------------------
mkdir -p "$dst/LOGS"

backup_and_copy() {
    local from="$1" to="$2"
    if [[ -f "$to" ]]; then
        cp -p "$to" "$to.backup_$stamp"
        echo "backup: $to -> $to.backup_$stamp"
    fi
    cp -p "$from" "$to"
    echo "copied: $from -> $to"
}

for f in "$src"/inlist*; do
    [[ -f "$f" ]] || continue
    backup_and_copy "$f" "$dst/$(basename "$f")"
done
backup_and_copy "$hist_src" "$dst/LOGS/history.data"

echo
echo "done"
