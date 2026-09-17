#pragma once

#include "coefficients.hpp"
#include "util.hpp"

namespace point_to_ellipse_series {

// Convert rational coefficient to type T.
template<typename T> T to(const mpq_class& c);

// Raise base to integer power.
template<typename T> T ipow(const T& base, int exp);

// sin(2*n * psi) — for multiple-angle series.
template<typename T> T isin_mul(const T& psi_v, int n);

// cos(2*n * psi) — for multiple-angle series.
template<typename T> T icos_mul(const T& psi_v, int n);

}  // namespace point_to_ellipse_series
