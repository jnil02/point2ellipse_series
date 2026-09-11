#pragma once

#include <mpreal.h>

#include "coefficients.hpp"
#include "util.hpp"

namespace point_to_ellipse_series {

// Convert rational coefficient to type T.
template<typename T>
T series_coeff(const mpq_class& c);

template<>
inline mpfr::mpreal series_coeff<mpfr::mpreal>(const mpq_class& c) {
	mpfr::mpreal num, den;
	mpfr_set_z(num.mpfr_ptr(), c.get_num().get_mpz_t(), mpfr::mpreal::get_default_rnd());
	mpfr_set_z(den.mpfr_ptr(), c.get_den().get_mpz_t(), mpfr::mpreal::get_default_rnd());
	return num / den;
}
template<> inline double series_coeff<double>(const mpq_class& c) {
	return c.get_d();
}

// Raise base to integer power.
template<typename T>
T series_pow(const T& base, int exp);

template<>
inline mpfr::mpreal series_pow<mpfr::mpreal>(const mpfr::mpreal& base, int exp) {
	return mpfr::pow(base, static_cast<long>(exp));
}
// FIXME(JO) The different versions should be put in different headers such that we don't need headers for both types if only one is needed.
template<> inline double series_pow<double>(const double& b, int e) {
	return std::pow(b, e);
}

// sin(2*n * psi) — for multiple-angle series.
template<typename T>
T series_sin_mul(const T& psi_v, int n);

template<>
inline mpfr::mpreal series_sin_mul<mpfr::mpreal>(const mpfr::mpreal& psi_v, int n) {
	return mpfr::sin(mpfr::mpreal(2 * n) * psi_v);
}

// cos(2*n * psi) — for multiple-angle series.
template<typename T>
T series_cos_mul(const T& psi_v, int n);

template<>
inline mpfr::mpreal series_cos_mul<mpfr::mpreal>(const mpfr::mpreal& psi_v, int n) {
	return mpfr::cos(mpfr::mpreal(2 * n) * psi_v);
}

}  // namespace point_to_ellipse_series
