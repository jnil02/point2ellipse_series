#pragma once
#include "coefficients_evo.hpp"

namespace point_to_ellipse_series::detail {

// Intermediate quantities and not part of coefficient-API,
// Kept accessible for the micro-benchmarks.
mpq_class a_mr(int m, int r);
mpq_class B_rt(int r, int t);
mpq_class C_mt(int m, int t);
mpq_class R(int n, int k, int l, int i);

} // namespace point_to_ellipse_series::detail