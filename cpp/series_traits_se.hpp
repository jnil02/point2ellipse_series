#pragma once

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/pow.h>
#include <symengine/rational.h>

#include "series_traits.hpp"

namespace point_to_ellipse_series {

inline SymEngine::Expression mpq_to_expr(const mpq_class& d) {
	SymEngine::rational_class q(d.get_mpq_t());
	return SymEngine::Expression(SymEngine::Rational::from_mpq(std::move(q)));
}

template<>
inline SymEngine::Expression series_coeff<SymEngine::Expression>(const mpq_class& c) {
	return mpq_to_expr(c);
}

template<>
inline SymEngine::Expression series_pow<SymEngine::Expression>(const SymEngine::Expression& base, int exp) {
	return SymEngine::pow(base, SymEngine::Expression(exp));
}

template<>
inline SymEngine::Expression series_sin_mul<SymEngine::Expression>(const SymEngine::Expression& psi_v, int n) {
	return SymEngine::sin(SymEngine::Expression(2 * n) * psi_v);
}

template<>
inline SymEngine::Expression series_cos_mul<SymEngine::Expression>(const SymEngine::Expression& psi_v, int n) {
	return SymEngine::cos(SymEngine::Expression(2 * n) * psi_v);
}

}  // namespace point_to_ellipse_series
