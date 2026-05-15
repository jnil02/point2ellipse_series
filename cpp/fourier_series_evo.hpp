#pragma once

/*
 * Symbolic inside-evolute series expansions for the point-to-ellipse relation.
 */

#include "fourier_series.hpp"
#include "coefficients_evo.hpp"

using point_to_ellipse_series::d_phi_evo;
using point_to_ellipse_series::c_phi_evo;
using point_to_ellipse_series::c_sin_phi_evo;
using point_to_ellipse_series::c_cos_phi_evo;
using point_to_ellipse_series::ch_evo;
using point_to_ellipse_series::dh_evo;

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
