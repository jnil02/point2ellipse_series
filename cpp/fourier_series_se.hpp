#pragma once

#include <symengine/expression.h>
#include <symengine/ntheory.h>

#include "fourier_series.hpp"
#include "series_traits_se.hpp"
#include "symbols.hpp"

using SymEngine::Expression;
using SymEngine::pow;
using SymEngine::sin;
using SymEngine::cos;

using point_to_ellipse_series::powm1;

inline Expression sigma(int J, const Expression& delta) {
	Expression d(0);
	for (int j = 0; j < J; j += 2) {
		int j2 = j / 2;
		d += Expression(SymEngine::div(SymEngine::integer(point_to_ellipse_series::powm1(j2)),
									   SymEngine::factorial(j)))
			 * pow(delta, Expression(j2));
	}
	return d;
}

inline Expression tau(int J, const Expression& omega, const Expression& delta) {
	Expression d(0);
	for (int j = 1; j < J; j += 2) {
		int j2 = (j - 1) / 2;
		d += Expression(SymEngine::div(SymEngine::integer(point_to_ellipse_series::powm1(j2)),
									   SymEngine::factorial(j)))
			 * pow(delta, Expression(j2));
	}
	return omega * d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_in_sin_pow(int N, int K) {
	return phi_in_sin_pow<Expression>(N, K, sin_psi, varrho, e2);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_in_sin_pow2(int N, int K) {
	return phi_in_sin_pow2<Expression>(N, K, sin_psi, varrho, e2);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_in_sin_mul(int N, int K, int L) {
	return phi_in_sin_mul<Expression>(N, K, L, psi, varrho, e2);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression sin_phi_in_sin_pow(int N, int K) {
	return sin_phi_in_sin_pow<Expression>(N, K, sin_psi, varrho, e2);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression sin_phi_in_cos_mul(int N, int K, int L) {
	return sin_phi_in_cos_mul<Expression>(N, K, L, psi, varrho, e2);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression cos_phi_in_sin_pow(int N, int K) {
	return cos_phi_in_sin_pow<Expression>(N, K, sin_psi, varrho, e2);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression cos_phi_in_cos_mul(int N, int K, int L) {
	return cos_phi_in_cos_mul<Expression>(N, K, L, psi, varrho, e2);
}

/** Series expansion of sin(phi) / sin(psi) - 1 in sin-powers but with common sigma/tau components. */
inline Expression sin_phi_in_sin_pow2(int N, int K, int J) {
	Expression sin_psi2 = pow(sin_psi, 2);
	Expression o = phi_in_sin_pow(N, K);
	Expression d = sin_psi2 * (Expression(1) - sin_psi2) * pow(o, 2);
	Expression s = sigma(J, d);
	Expression t = tau(J, o, d);
	return s + (Expression(1) - sin_psi2) * t;
}

/** Series expansion of cos(phi) / cos(psi) - 1 in sin-powers but with common sigma/tau components. */
inline Expression cos_phi_in_sin_pow2(int N, int K, int J) {
	Expression sin_psi2 = pow(sin_psi, 2);
	Expression o = phi_in_sin_pow(N, K);
	Expression d = sin_psi2 * (Expression(1) - sin_psi2) * pow(o, 2);
	Expression s = sigma(J, d);
	Expression t = tau(J, o, d);
	return s - sin_psi2 * t;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_in_sin_pow(int N, int K) {
	return h_in_sin_pow<Expression>(N, K, sin_psi, varrho, e2);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_in_cos_mul(int N, int K, int L) {
	return h_in_cos_mul<Expression>(N, K, L, psi, varrho, e2);
}
