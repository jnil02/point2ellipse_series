#pragma once

#include <symengine/expression.h>
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
	return mpfr::mpreal(c.num) / mpfr::mpreal(c.den);
}

// Raise base to integer power.
inline SymEngine::Expression series_pow(const SymEngine::Expression& base, int exp) {
	return SymEngine::pow(base, SymEngine::Expression(exp));
}

inline mpfr::mpreal series_pow(const mpfr::mpreal& base, int exp) {
	return mpfr::pow(base, static_cast<long>(exp));
}

}  // namespace point_to_ellipse_series
