#!/usr/bin/env bash
# Run cold-call benchmarks for key coefficient functions at increasing orders.
# Output is tab-separated: function  args...  wall_ms
#
# Usage:
#   ./run_benchmarks.sh [path/to/bench_coefficients]
#
# If no path is given, looks for the binary next to this script's build
# output (../build/bench_coefficients relative to this script).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH="${1:-"$SCRIPT_DIR/../build/bench_coefficients"}"

if [[ ! -x "$BENCH" ]]; then
    echo "bench_coefficients binary not found at: $BENCH" >&2
    echo "Build it first: cd cpp/build && make bench_coefficients" >&2
    exit 1
fi

echo -e "function\targs\twall_ms"
echo -e "--------\t----\t-------"

# a_mr(m, r) — moderate growth, good warm-up
for m in 2 4 6 8 10; do
    r=$((m / 2))
    "$BENCH" a_mr $m $r
done

echo ""

# d_phi_evo(n, k, l) — shallowest coefficient.
for n in 2 4 6 8 10; do
    k=$n
    l=$n
    "$BENCH" d_phi_evo $n $k $l
done

echo ""

# d_h_evo(n, k, l) — deepest coefficient.
for n in 2 4 6 8; do
    k=$n
    l=$n
    "$BENCH" d_h_evo $n $n $n
done
