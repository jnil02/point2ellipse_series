#pragma once

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/pow.h>
#include <mpreal.h>

#include "coefficients.hpp"
#include "util.hpp"

namespace point_to_ellipse_series {

// Convert rational coefficient rc to type T.
template<typename T>
T series_coeff(const rc& c);

template<>
inline SymEngine::Expression series_coeff<SymEngine::Expression>(const rc& c) {
	return rc_expr(c);
}

template<>
inline mpfr::mpreal series_coeff<mpfr::mpreal>(const rc& c) {
	mpfr::mpreal num, den;
	mpfr_set_z(num.mpfr_ptr(), c.num.get_mpz_t(), mpfr::mpreal::get_default_rnd());
	mpfr_set_z(den.mpfr_ptr(), c.den.get_mpz_t(), mpfr::mpreal::get_default_rnd());
	return num / den;
}

// Raise base to integer power.
inline SymEngine::Expression series_pow(const SymEngine::Expression& base, int exp) {
	return SymEngine::pow(base, SymEngine::Expression(exp));
}

inline mpfr::mpreal series_pow(const mpfr::mpreal& base, int exp) {
	return mpfr::pow(base, static_cast<long>(exp));
}

// sin(2*n * psi) — for multiple-angle series.
inline SymEngine::Expression series_sin_mul(const SymEngine::Expression& psi_v, int n) {
	return SymEngine::sin(SymEngine::Expression(2 * n) * psi_v);
}

inline mpfr::mpreal series_sin_mul(const mpfr::mpreal& psi_v, int n) {
	return mpfr::sin(mpfr::mpreal(2 * n) * psi_v);
}

// cos(2*n * psi) — for multiple-angle series.
inline SymEngine::Expression series_cos_mul(const SymEngine::Expression& psi_v, int n) {
	return SymEngine::cos(SymEngine::Expression(2 * n) * psi_v);
}

inline mpfr::mpreal series_cos_mul(const mpfr::mpreal& psi_v, int n) {
	return mpfr::cos(mpfr::mpreal(2 * n) * psi_v);
}

}  // namespace point_to_ellipse_series
