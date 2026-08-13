#!/bin/bash
# Build-time smoke tests for the monolithic image -- fail the build loudly
# rather than pushing something broken to quay.io. Run from container/Dockerfile.
#
# Every per-rule tool is checked INSIDE whichever pre-built conda env provides
# it, not on the base PATH: none of them are on the base PATH, and that's the
# point of the per-rule env layout. Envs are located by searching the prefix
# rather than by hardcoded directory names, because Snakemake names each env
# by a content hash that this script has no business predicting.

set -euo pipefail

PREFIX="$(python -c 'import phables, os; print(os.path.join(os.path.dirname(phables.__file__), "workflow", "conda"))')"
echo "conda prefix: $PREFIX"
test -d "$PREFIX"

echo "=== phables CLI ==="
phables --version
phables run -h > /dev/null
phables install -h > /dev/null
echo "CLI OK"

echo "=== the --container / --prostt5-container flags must be GONE ==="
if phables run -h 2>&1 | grep -qE '\-\-(prostt5-)?container'; then
    echo "ERROR: a container flag is still present in the CLI" >&2
    exit 1
fi
echo "confirmed absent"

echo "=== pre-built env count ==="
n=$(find "$PREFIX" -maxdepth 1 -mindepth 1 -type d | wc -l)
echo "found $n env directories"
find "$PREFIX" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;
# 6 distinct env files are reachable at minimum (coverm, genecall, smg, mmseqs,
# phables, curl) before counting foldseek/prostt5-rocm/prostt5-cpu/phylotree.
test "$n" -ge 6

# Finds an executable in any pre-built env; fails if no env provides it.
check_bin() {
    local want="$1" e
    for e in "$PREFIX"/*/; do
        if [ -x "${e}bin/${want}" ]; then
            echo "OK: $want -> $e"
            return 0
        fi
    done
    echo "MISSING from every pre-built env: $want" >&2
    return 1
}

echo "=== per-rule binaries ==="
check_bin minimap2
check_bin samtools
check_bin coverm
check_bin mmseqs
check_bin foldseek
check_bin hmmsearch
# the bioconda package's real binary name, confirmed against genes.smk's own
# invocation -- NOT `fraggenescan`
check_bin run_FragGeneScan.pl
check_bin mafft
check_bin curl

echo "=== torch must come from the prostt5-rocm env, not the base image ==="
# The base image ships its own system-python torch 2.7.1, which the workflow
# does NOT use. predict_3di runs inside the prostt5-rocm conda env, whose torch
# is 2.9.1+rocm6.3 (workflow/envs/prostt5-rocm.yaml). If nothing here reports
# 2.9.1+rocm, that env didn't build properly and predict_3di would be running
# on the wrong stack.
found=0
for e in "$PREFIX"/*/; do
    if [ -x "${e}bin/python" ]; then
        if "${e}bin/python" - <<'PY' 2>/dev/null
import sys
try:
    import torch
except Exception:
    sys.exit(1)
sys.exit(0 if torch.__version__.startswith("2.9.1") and "rocm" in torch.__version__ else 1)
PY
        then
            echo "OK: torch 2.9.1+rocm in $e"
            "${e}bin/python" -c "import pholdlib; print('pholdlib OK')"
            found=1
            break
        fi
    fi
done
if [ "$found" -ne 1 ]; then
    echo "ERROR: no pre-built env provides torch 2.9.1+rocm" >&2
    exit 1
fi

echo "=== all image tests passed ==="
