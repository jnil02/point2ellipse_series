#pragma once

/*
 * Symbolic series expansions for the point-to-ellipse relation.
 */

#include <symengine/expression.h>

#include "coefficients.hpp"
#include "polynomials.hpp"
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

using point_to_ellipse_series::d_phi_evo;
using point_to_ellipse_series::c_phi_evo;
using point_to_ellipse_series::c_sin_phi_evo;
using point_to_ellipse_series::c_cos_phi_evo;
using point_to_ellipse_series::ch_evo;
using point_to_ellipse_series::dh_evo;

using point_to_ellipse_series::rc_expr;
using point_to_ellipse_series::sigma;
using point_to_ellipse_series::tau;

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
				d += rc_expr(d_phi(n, k, l)) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n);
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
				d += rc_expr(d_phi2(n, k, l)) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n + 1);
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
				d += rc_expr(c_phi(n, k, l)) * pow(e2, l) * pow(varrho, k) * sin(Expression(2 * n) * psi);
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
				d += rc_expr(d_sin(n, k, l)) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n);
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
				d += rc_expr(c_sin(n, k, l)) * pow(e2, l) * pow(varrho, k) * cos(Expression(2 * n) * psi);
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
				d += rc_expr(d_cos(n, k, l)) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n);
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
				d += rc_expr(c_cos(n, k, l)) * pow(e2, l) * pow(varrho, k) * cos(Expression(2 * n) * psi);
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
				d += rc_expr(d_h(n, k, l)) * pow(e2, l) * pow(varrho, k) * pow(sin_psi, 2 * n);
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
				d += rc_expr(c_h(n, k, l)) * pow(e2, l) * pow(varrho, k) * cos(Expression(2 * n) * psi);
	return d;
}

// ---------------------------------------------------------------------------
// Inside-evolute series
// ---------------------------------------------------------------------------

/** Inside-evolute series for (phi - sgn*pi/2) / (sgn*|cos(psi)|) in sin-powers (dense).
 *
 * @param N sin power series truncation order.
 * @param K rho_ae2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression phi_evo_sin_pow_dense(int N, int K) {
	Expression d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 0; k <= K; ++k)
			for (int l = 0; l <= n / 2 + k; ++l)
				d += rc_expr(d_phi_evo(n, k, l))
					 * pow(sin_psi, n)
					 * pow(rho_ae2, n + 1 + 2 * k)
					 * pow(b_a, (n % 2) + 1 + 2 * l);
	return d;
}

/** Inside-evolute series for (phi - sgn*pi/2) / (sgn*|cos(psi)|) in sin-powers (sparse).
 *
 * @param N sin power series truncation order.
 * @param K rho_ae2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression phi_evo_sin_pow(int N, int K) {
	Expression d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = n + 1; k <= K; ++k)
			for (int l = 1; l <= k; ++l)
				d += rc_expr(c_phi_evo(n, k, l))
					 * cos_psi * pow(sin_psi, n)
					 * pow(rho_ae2, k) * pow(b_a, l);
	return d;
}

/** Inside-evolute series for sin(phi) in sin-powers (dense).
 *
 * @param N sin power series truncation order.
 * @param K rho_ae2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression sin_phi_evo_dense(int N, int K) {
	Expression s(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 0; n + 2 * k <= K; ++k)
			for (int l = 1; l <= n / 2 + k; ++l)
				s += rc_expr(c_sin_phi_evo(n, n + 2 * k, 2 * l + (n % 2)))
					 * pow(b_a, 2 * l + (n % 2))
					 * pow(rho_ae2, n + 2 * k)
					 * pow(sin_psi, n);
	return s;
}

/** Inside-evolute series for cos(phi) / |cos(psi)| in sin-powers (dense).
 *
 * @param N sin power series truncation order.
 * @param K rho_ae2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression cos_phi_evo_dense(int N, int K) {
	Expression s(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 0; n + 1 + 2 * k <= K; ++k)
			for (int l = n % 2; l <= (n + 1) / 2 + k; ++l)
				s += rc_expr(c_cos_phi_evo(n, n + 1 + 2 * k, 2 * l + 1 - (n % 2)))
					 * pow(b_a, 2 * l + 1 - (n % 2))
					 * pow(rho_ae2, n + 1 + 2 * k)
					 * pow(sin_psi, n);
	return s;
}

/** Inside-evolute series for h/a in sin-powers.
 *
 * @param N sin power series truncation order.
 * @param K rho_ae2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression h_a_evo(int N, int K) {
	Expression s(0);
	for (int n = 0; n <= N; ++n)
		for (int k = n; k <= K; ++k)
			for (int l = 0; l <= k + 1; ++l)
				s += rc_expr(ch_evo(n, k, l))
					 * pow(b_a, l) * pow(rho_ae2, k) * pow(sin_psi, n);
	return s;
}

/** Inside-evolute series for h/a in sin-powers (dense).
 *
 * @param N sin power series truncation order.
 * @param K rho_ae2 power series truncation order.
 * @return Symbolic expression for the truncated series.
 */
inline Expression h_a_evo_dense(int N, int K) {
	Expression s(0);
	for (int n = 0; n <= N; ++n) {
		int sn = n % 2;
		for (int k = 0; n + 2 * k <= K; ++k)
			for (int l = 0; l <= k + (n + 1) / 2; ++l)
				s += rc_expr(dh_evo(n, k, l))
					 * pow(b_a, 1 - sn + 2 * l)
					 * pow(rho_ae2, n + 2 * k)
					 * pow(sin_psi, n);
	}
	return s;
}