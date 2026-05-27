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
using point_to_ellipse_series::c_h_evo;
using point_to_ellipse_series::d_h_evo;

/** Inside-evolute series for (phi - sgn*pi/2) / (sgn*|cos(psi)|) in sin-powers (dense).
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho_ae2 power series truncation order.
 * @param sin_psi_v   Value/expression for |sin(psi)|.
 * @param rho_ae2_v   Value/expression for rho/(a*e²).
 * @param b_a_v       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T phi_evo_sin_pow_dense(int N, int K,
							   const T& sin_psi_v, const T& rho_ae2_v, const T& b_a_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 0; k <= K; ++k)
			for (int l = 0; l <= n / 2 + k; ++l)
				d = d + series_coeff<T>(d_phi_evo(n, k, l))
						* series_pow(sin_psi_v, n)
						* series_pow(rho_ae2_v, n + 1 + 2 * k)
						* series_pow(b_a_v, (n % 2) + 1 + 2 * l);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_evo_sin_pow_dense(int N, int K) {
	return phi_evo_sin_pow_dense<Expression>(N, K, sin_psi, rho_ae2, b_a);
}

template<typename T>
inline T phi_evo_sin_pow_dense_m(int M,
								 const T& sin_psi_v, const T& rho_ae2_v, const T& b_a_v) {
	T d(0);
	for (int m = 1; m <= M; ++m) {
		const T   rho_m  = series_pow(rho_ae2_v, m);
		const int L      = (m - 1) / 2;   // floor((m-1)/2): upper bound for k and l
		const int parity = (m - 1) % 2;   // n%2 == (m-1)%2 for all valid n at this m
		for (int k = 0; k <= L; ++k) {
			const int n     = m - 1 - 2 * k;
			const T   sin_n = series_pow(sin_psi_v, n);
			for (int l = 0; l <= L; ++l)
				d = d + series_coeff<T>(d_phi_evo(n, k, l))
						* sin_n
						* rho_m
						* series_pow(b_a_v, parity + 1 + 2 * l);
		}
	}
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_evo_sin_pow_dense_m(int K) {
	return phi_evo_sin_pow_dense_m<Expression>(K, sin_psi, rho_ae2, b_a);
}


/** Inside-evolute series for (phi - sgn*pi/2) / (sgn*|cos(psi)|) in sin-powers (sparse).
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho_ae2 power series truncation order.
 * @param sin_psi_v   Value/expression for |sin(psi)|.
 * @param cos_psi_v   Value/expression for |cos(psi)|.
 * @param rho_ae2_v   Value/expression for rho/(a*e²).
 * @param b_a_v       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T phi_evo_sin_pow(int N, int K,
						 const T& sin_psi_v, const T& cos_psi_v,
						 const T& rho_ae2_v, const T& b_a_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = n + 1; k <= K; ++k)
			for (int l = 1; l <= k; ++l)
				d = d + series_coeff<T>(c_phi_evo(n, k, l))
						* cos_psi_v
						* series_pow(sin_psi_v, n)
						* series_pow(rho_ae2_v, k)
						* series_pow(b_a_v, l);
	return d;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_evo_sin_pow(int N, int K) {
	return phi_evo_sin_pow<Expression>(N, K, sin_psi, cos_psi, rho_ae2, b_a);
}

/** Inside-evolute series for sin(phi) in sin-powers (dense).
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho_ae2 power series truncation order.
 * @param sin_psi_v   Value/expression for |sin(psi)|.
 * @param rho_ae2_v   Value/expression for rho/(a*e²).
 * @param b_a_v       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T sin_phi_evo_dense(int N, int K,
						   const T& sin_psi_v, const T& rho_ae2_v, const T& b_a_v) {
	T s(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 0; n + 2 * k <= K; ++k)
			for (int l = 1; l <= n / 2 + k; ++l)
				s = s + series_coeff<T>(c_sin_phi_evo(n, n + 2 * k, 2 * l + (n % 2)))
						* series_pow(b_a_v, 2 * l + (n % 2))
						* series_pow(rho_ae2_v, n + 2 * k)
						* series_pow(sin_psi_v, n);
	return s;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression sin_phi_evo_dense(int N, int K) {
	return sin_phi_evo_dense<Expression>(N, K, sin_psi, rho_ae2, b_a);
}

/** Inside-evolute series for cos(phi) / |cos(psi)| in sin-powers (dense).
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho_ae2 power series truncation order.
 * @param sin_psi_v   Value/expression for |sin(psi)|.
 * @param rho_ae2_v   Value/expression for rho/(a*e²).
 * @param b_a_v       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T cos_phi_evo_dense(int N, int K,
						   const T& sin_psi_v, const T& rho_ae2_v, const T& b_a_v) {
	T s(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 0; n + 1 + 2 * k <= K; ++k)
			for (int l = n % 2; l <= (n + 1) / 2 + k; ++l)
				s = s + series_coeff<T>(c_cos_phi_evo(n, n + 1 + 2 * k, 2 * l + 1 - (n % 2)))
						* series_pow(b_a_v, 2 * l + 1 - (n % 2))
						* series_pow(rho_ae2_v, n + 1 + 2 * k)
						* series_pow(sin_psi_v, n);
	return s;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression cos_phi_evo_dense(int N, int K) {
	return cos_phi_evo_dense<Expression>(N, K, sin_psi, rho_ae2, b_a);
}

/** Inside-evolute series for h/a in sin-powers.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho_ae2 power series truncation order.
 * @param sin_psi_v   Value/expression for |sin(psi)|.
 * @param rho_ae2_v   Value/expression for rho/(a*e²).
 * @param b_a_v       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T h_a_evo(int N, int K,
				 const T& sin_psi_v, const T& rho_ae2_v, const T& b_a_v) {
	T s(0);
	for (int n = 0; n <= N; ++n)
		for (int k = n; k <= K; ++k)
			for (int l = 0; l <= k + 1; ++l)
				s = s + series_coeff<T>(c_h_evo(n, k, l))
						* series_pow(b_a_v, l)
						* series_pow(rho_ae2_v, k)
						* series_pow(sin_psi_v, n);
	return s;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_a_evo(int N, int K) {
	return h_a_evo<Expression>(N, K, sin_psi, rho_ae2, b_a);
}

/** Inside-evolute series for h/a in sin-powers (dense).
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho_ae2 power series truncation order.
 * @param sin_psi_v   Value/expression for |sin(psi)|.
 * @param rho_ae2_v   Value/expression for rho/(a*e²).
 * @param b_a_v       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T h_a_evo_dense(int N, int K,
					   const T& sin_psi_v, const T& rho_ae2_v, const T& b_a_v) {
	T s(0);
	for (int n = 0; n <= N; ++n) {
		int sn = n % 2;
		for (int k = 0; n + 2 * k <= K; ++k)
			for (int l = 0; l <= k + (n + 1) / 2; ++l)
				s = s + series_coeff<T>(d_h_evo(n, k, l))
						* series_pow(b_a_v, 1 - sn + 2 * l)
						* series_pow(rho_ae2_v, n + 2 * k)
						* series_pow(sin_psi_v, n);
	}
	return s;
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_a_evo_dense(int N, int K) {
	return h_a_evo_dense<Expression>(N, K, sin_psi, rho_ae2, b_a);
}
