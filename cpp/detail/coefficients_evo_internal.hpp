#pragma once

#include <gmpxx.h>

namespace point_to_ellipse_series::detail {

// Intermediate quantities, not part of the coefficient API.
// Kept accessible for the micro-benchmarks.

/// Helper coefficient for the ε²·N/a reduction.
mpq_class a_mr(int m, int r);
/// Helper coefficient for the ε²·N/a reduction.
mpq_class B_rt(int r, int t);
/// Helper coefficient for the ε²·N/a reduction.
mpq_class C_mt(int m, int t);
/// Intermediate alternating sum feeding c_N_evo.
mpq_class R(int k, int l, int n, int i);

} // namespace point_to_ellipse_series::detail
