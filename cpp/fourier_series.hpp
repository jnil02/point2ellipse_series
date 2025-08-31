#pragma once

/*
 * Symbolic series expansions for the point-to-ellipse relation.
 */

#include <symengine/expression.h>

#include "coefficients.hpp"
#include "symbols.hpp"

using SymEngine::Expression;
using SymEngine::pow;
using SymEngine::sin;
using SymEngine::cos;

/** Series expansion of (phi - psi) / (sin(psi) * cos(psi)) in sin-powers.
 *
 * @param N sin power series truncation order.
 * @param K rho power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression phi_in_sin_pow(int N, int K) {
	Expression d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(k, n + 1); l <= k + n; ++l)
				d += d_phi(n, k, l) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n);
	return d;
}

/** Series expansion of (phi - psi) in sin-powers.
 *
 * Note, this series has poor convergence and is only implemented to demonstrate
 * this. It should not be used in practice.
 *
 * @param N sin power series truncation order.
 * @param K rho power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression phi_in_sin_pow2(int N, int K) {
	Expression d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = k; l <= k + n; ++l)
				d += d_phi2(n, k, l) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n + 1);
	return d;
}

/** Series expansion of (phi - psi) in sin multiples.
 *
 * @param N sin multiple series truncation order.
 * @param K rho power series truncation order.
 * @param L e2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
Expression phi_in_sin_mul(int N, int K, int L) {
	Expression d(0);
	for (int n = 1; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(n, k); l <= L; ++l)
				d += c_phi(n, k, l) * pow(e2, l) * pow(varrho, k) * sin(Expression(2 * n) * psi);
	return d;
}

/** Series expansion of sin(phi) / sin(psi) - 1 in sin-powers.
 *
 * @param N sin power series truncation order.
 * @param K rho power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression sin_phi_in_sin_pow(int N, int K) {
	Expression d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(k, n); l <= n + k; ++l)
				d += d_sin(n, k, l) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n);
	return d;
}

/** Series expansion of sin(phi) / sin(psi) - 1 in cos multiples.
 *
 * @param N sin multiple series truncation order.
 * @param K rho power series truncation order.
 * @param L e2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression sin_phi_in_cos_mul(int N, int K, int L) {
	Expression d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(n, k); l <= L; ++l)
				d += c_sin(n, k, l) * pow(e2, l) * pow(varrho, k) * cos(Expression(2 * n) * psi);
	return d;
}

/** Series expansion of cos(phi) / cos(psi) - 1 in sin-powers.
 *
 * @param N sin power series truncation order.
 * @param K rho power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression cos_phi_in_sin_pow(int N, int K) {
	Expression d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(k, n); l < n + k; ++l)
				d += d_cos(n, k, l) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n);
	return d;
}

/** Series expansion of cos(phi) / cos(psi) - 1 in cos multiples.
 *
 * @param N sin multiple series truncation order.
 * @param K rho power series truncation order.
 * @param L e2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression cos_phi_in_cos_mul(int N, int K, int L) {
	Expression d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(n, k); l <= L; ++l)
				d += c_cos(n, k, l) * pow(e2, l) * pow(varrho, k) * cos(Expression(2 * n) * psi);
	return d;
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
 * @param N sin power series truncation order.
 * @param K rho power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression h_in_sin_pow(int N, int K) {
	Expression d(0);
	for (int n = 1; n <= N; ++n)
		for (int k = 0; k <= K; ++k)
			for (int l = std::max(k + 1, n); l <= n + k; ++l)
				d += d_h(n, k, l) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n);
	return d;
}

/** Series expansion of (h + a + rho) / a in cos multiples.
 *
 * @param N sin multiples series truncation order.
 * @param K rho power series truncation order.
 * @param L e2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression h_in_cos_mul(int N, int K, int L) {
	Expression d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 0; k <= K; ++k)
			for (int l = std::max(n, k + 1); l <= L; ++l)
				d += c_h(n, k, l) * pow(e2, l) * pow(varrho, k) * cos(Expression(2 * n) * psi);
	return d;
}
