#!/usr/bin/env bash
# Run cold-call benchmarks for coefficient functions and series evaluation.
# Output is tab-separated: function  args...  wall_ms
#
# Usage:
#   ./run_benchmarks.sh [build_dir]
#
# If no build_dir is given, defaults to ../build relative to this script.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${1:-"$SCRIPT_DIR/../build"}"
BENCH_COEFF="$BUILD/bench_coefficients"
BENCH_SERIES="$BUILD/bench_series_evo"

for bin in "$BENCH_COEFF" "$BENCH_SERIES"; do
    if [[ ! -x "$bin" ]]; then
        echo "Binary not found: $bin" >&2
        echo "Build with: cd cpp/build && make bench_coefficients bench_series_evo" >&2
        exit 1
    fi
done

# Representative point: ellipse with b/a=0.97 (~WGS-84 flattening), moderate rho.
SIN_PSI=0.5
RHO_AE2=0.3
B_A=0.97

echo "=== Coefficient benchmarks (cold, one call per process) ==="
echo -e "function\targs\twall_ms"

echo "--- a_mr ---"
for m in 2 4 6 8 10; do
    "$BENCH_COEFF" a_mr $m $((m/2))
done

echo "--- d_phi_evo ---"
for n in 2 4 6 8 10; do
    "$BENCH_COEFF" d_phi_evo $n $n $n
done

echo "--- d_h_evo ---"
for n in 2 4 6; do
    "$BENCH_COEFF" d_h_evo $n $n $n
done

echo ""
echo "=== Series evaluation benchmarks (mpreal, includes warm coefficient cache) ==="
echo -e "series\tN\tK\twall_ms"

echo "--- phi_evo_dense_m (the sweep series) ---"
for M in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
    "$BENCH_SERIES" phi_evo_dense_m $M $M $SIN_PSI $RHO_AE2 $B_A
done

echo "--- phi_evo_dense ---"
for N in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
    "$BENCH_SERIES" phi_evo_dense $N $N $SIN_PSI $RHO_AE2 $B_A
done

echo "--- sin_phi_evo_dense ---"
for N in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
    "$BENCH_SERIES" sin_phi_evo_dense $N $N $SIN_PSI $RHO_AE2 $B_A
done

echo "--- cos_phi_evo_dense ---"
for N in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
    "$BENCH_SERIES" cos_phi_evo_dense $N $N $SIN_PSI $RHO_AE2 $B_A
done

echo "--- h_a_evo_dense ---"
for N in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
    "$BENCH_SERIES" h_a_evo_dense $N $N $SIN_PSI $RHO_AE2 $B_A
done
