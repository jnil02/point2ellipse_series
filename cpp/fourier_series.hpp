#pragma once

/*
 * Symbolic series expansions for the point-to-ellipse relation.
 */

#include <symengine/expression.h>
#include <symengine/ntheory.h>

#include "coefficients.hpp"
#include "polynomials.hpp"
#include "series_traits.hpp"
#include "symbols.hpp"

using SymEngine::Expression;
using SymEngine::pow;
using SymEngine::sin;
using SymEngine::cos;

using point_to_ellipse_series::d_phi;
using point_to_ellipse_series::d_phi2;
using point_to_ellipse_series::c_phi;
using point_to_ellipse_series::d_sin;
using point_to_ellipse_series::c_sin;
using point_to_ellipse_series::d_cos;
using point_to_ellipse_series::c_cos;
using point_to_ellipse_series::d_h;
using point_to_ellipse_series::c_h;

using point_to_ellipse_series::series_coeff;
using point_to_ellipse_series::series_pow;
using point_to_ellipse_series::series_sin_mul;
using point_to_ellipse_series::series_cos_mul;

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

/** Series expansion of (phi - psi) / (sin(psi) * cos(psi)) in sin-powers.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for rho/a.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T phi_in_sin_pow(int N, int K,
						const T& sin_psi_v, const T& varrho_v, const T& e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(k, n + 1); l <= k + n; ++l)
				d = d + series_coeff<T>(d_phi(n, k, l))
						* series_pow(e2_v, l)
						* series_pow(varrho_v, k)
						* series_pow(sin_psi_v, 2 * n);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_in_sin_pow(int N, int K) {
	return phi_in_sin_pow<Expression>(N, K, sin_psi, varrho, e2);
}

/** Series expansion of (phi - psi) in sin-powers.
 *
 * Note, this series has poor convergence and is only implemented to demonstrate
 * this. It should not be used in practice.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for rho/a.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T phi_in_sin_pow2(int N, int K,
						 const T& sin_psi_v, const T& varrho_v, const T& e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = k; l <= k + n; ++l)
				d = d + series_coeff<T>(d_phi2(n, k, l))
						* series_pow(e2_v, l)
						* series_pow(varrho_v, k)
						* series_pow(sin_psi_v, 2 * n + 1);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_in_sin_pow2(int N, int K) {
	return phi_in_sin_pow2<Expression>(N, K, sin_psi, varrho, e2);
}

/** Series expansion of (phi - psi) in sin multiples.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin multiple series truncation order.
 * @param K     rho power series truncation order.
 * @param L     e2 power series truncation order.
 * @param psi_v       Value/expression for psi.
 * @param varrho_v    Value/expression for rho/a.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T phi_in_sin_mul(int N, int K, int L,
						const T& psi_v, const T& varrho_v, const T& e2_v) {
	T d(0);
	for (int n = 1; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(n, k); l <= L; ++l)
				d = d + series_coeff<T>(c_phi(n, k, l))
						* series_pow(e2_v, l)
						* series_pow(varrho_v, k)
						* series_sin_mul(psi_v, n);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_in_sin_mul(int N, int K, int L) {
	return phi_in_sin_mul<Expression>(N, K, L, psi, varrho, e2);
}

/** Series expansion of sin(phi) / sin(psi) - 1 in sin-powers.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for rho/a.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T sin_phi_in_sin_pow(int N, int K,
							const T& sin_psi_v, const T& varrho_v, const T& e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(k, n); l <= n + k; ++l)
				d = d + series_coeff<T>(d_sin(n, k, l))
						* series_pow(e2_v, l)
						* series_pow(varrho_v, k)
						* series_pow(sin_psi_v, 2 * n);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression sin_phi_in_sin_pow(int N, int K) {
	return sin_phi_in_sin_pow<Expression>(N, K, sin_psi, varrho, e2);
}

/** Series expansion of sin(phi) / sin(psi) - 1 in cos multiples.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin multiple series truncation order.
 * @param K     rho power series truncation order.
 * @param L     e2 power series truncation order.
 * @param psi_v       Value/expression for psi.
 * @param varrho_v    Value/expression for rho/a.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T sin_phi_in_cos_mul(int N, int K, int L,
							const T& psi_v, const T& varrho_v, const T& e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(n, k); l <= L; ++l)
				d = d + series_coeff<T>(c_sin(n, k, l))
						* series_pow(e2_v, l)
						* series_pow(varrho_v, k)
						* series_cos_mul(psi_v, n);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression sin_phi_in_cos_mul(int N, int K, int L) {
	return sin_phi_in_cos_mul<Expression>(N, K, L, psi, varrho, e2);
}

/** Series expansion of cos(phi) / cos(psi) - 1 in sin-powers.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for rho/a.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T cos_phi_in_sin_pow(int N, int K,
							const T& sin_psi_v, const T& varrho_v, const T& e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(k, n); l < n + k; ++l)
				d = d + series_coeff<T>(d_cos(n, k, l))
						* series_pow(e2_v, l)
						* series_pow(varrho_v, k)
						* series_pow(sin_psi_v, 2 * n);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression cos_phi_in_sin_pow(int N, int K) {
	return cos_phi_in_sin_pow<Expression>(N, K, sin_psi, varrho, e2);
}

/** Series expansion of cos(phi) / cos(psi) - 1 in cos multiples.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin multiple series truncation order.
 * @param K     rho power series truncation order.
 * @param L     e2 power series truncation order.
 * @param psi_v       Value/expression for psi.
 * @param varrho_v    Value/expression for rho/a.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T cos_phi_in_cos_mul(int N, int K, int L,
							const T& psi_v, const T& varrho_v, const T& e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(n, k); l <= L; ++l)
				d = d + series_coeff<T>(c_cos(n, k, l))
						* series_pow(e2_v, l)
						* series_pow(varrho_v, k)
						* series_cos_mul(psi_v, n);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression cos_phi_in_cos_mul(int N, int K, int L) {
	return cos_phi_in_cos_mul<Expression>(N, K, L, psi, varrho, e2);
}

/** Series expansion of sin(phi) / sin(psi) - 1 in sin-powers but with common sigma/tau components.
 *
 * @param N sin power series truncation order.
 * @param K rho power series truncation order.
 * @param J sigma/tau power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression sin_phi_in_sin_pow2(int N, int K, int J) {
	Expression sin_psi2 = pow(sin_psi, 2);
	Expression o = phi_in_sin_pow(N, K);
	Expression d = sin_psi2 * (Expression(1) - sin_psi2) * pow(o, 2);
	Expression s = sigma(J, d);
	Expression t = tau(J, o, d);
	return s + (Expression(1) - sin_psi2) * t;
}

/** Series expansion of cos(phi) / cos(psi) - 1 in sin-powers but with common sigma/tau components.
 *
 * @param N sin power series truncation order.
 * @param K rho power series truncation order.
 * @param J sigma/tau power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression cos_phi_in_sin_pow2(int N, int K, int J) {
	Expression sin_psi2 = pow(sin_psi, 2);
	Expression o = phi_in_sin_pow(N, K);
	Expression d = sin_psi2 * (Expression(1) - sin_psi2) * pow(o, 2);
	Expression s = sigma(J, d);
	Expression t = tau(J, o, d);
	return s - sin_psi2 * t;
}

/** Series expansion of (h + a + rho) / a in sin-powers.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for rho/a.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T h_in_sin_pow(int N, int K,
					  const T& sin_psi_v, const T& varrho_v, const T& e2_v) {
	T d(0);
	for (int n = 1; n <= N; ++n)
		for (int k = 0; k <= K; ++k)
			for (int l = std::max(k + 1, n); l <= n + k; ++l)
				d = d + series_coeff<T>(d_h(n, k, l))
						* series_pow(e2_v, l)
						* series_pow(varrho_v, k)
						* series_pow(sin_psi_v, 2 * n);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_in_sin_pow(int N, int K) {
	return h_in_sin_pow<Expression>(N, K, sin_psi, varrho, e2);
}

/** Series expansion of (h + a + rho) / a in cos multiples.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin multiples series truncation order.
 * @param K     rho power series truncation order.
 * @param L     e2 power series truncation order.
 * @param psi_v       Value/expression for psi.
 * @param varrho_v    Value/expression for rho/a.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T h_in_cos_mul(int N, int K, int L,
					  const T& psi_v, const T& varrho_v, const T& e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 0; k <= K; ++k)
			for (int l = std::max(n, k + 1); l <= L; ++l)
				d = d + series_coeff<T>(c_h(n, k, l))
						* series_pow(e2_v, l)
						* series_pow(varrho_v, k)
						* series_cos_mul(psi_v, n);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_in_cos_mul(int N, int K, int L) {
	return h_in_cos_mul<Expression>(N, K, L, psi, varrho, e2);
}
