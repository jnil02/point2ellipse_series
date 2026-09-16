#pragma once

/*
 * Symbolic inside-evolute series expansions for the point-to-ellipse relation.
 */

#include "fourier_series.hpp"
#include "coefficients_evo.hpp"

using point_to_ellipse_series::d_phi_evo2;
using point_to_ellipse_series::d_sin_phi_evo2;
using point_to_ellipse_series::d_cos_phi_evo_m;
using point_to_ellipse_series::c_phi_evo;
using point_to_ellipse_series::c_sin_phi_evo;
using point_to_ellipse_series::c_cos_phi_evo;
using point_to_ellipse_series::c_h_evo;
using point_to_ellipse_series::c_h_evo2;
using point_to_ellipse_series::d_h_evo3;

/** Inside-evolute series for (phi - sgn*pi/2) / (sgn*|cos(psi)|) in sin-powers (sparse).
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param L     sin power series truncation order.
 * @param K     rho_ae2 power series truncation order.
 * @param abs_sin_psi   Value/expression for |sin(psi)|.
 * @param cos_psi_v   Value/expression for |cos(psi)|.
 * @param rho_ae2   Value/expression for rho/(a*e²).
 * @param b_a       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T phi_evo_sparse(int L, int K,
						const T& abs_sin_psi, const T& rho_ae2, const T& b_a) {
	T d(0);
	for (int l = 0; l <= L; ++l)
		for (int k = l + 1; k <= K; ++k)
			for (int n = 1; n <= k; ++n)
				d = d + point_to_ellipse_series::series_coeff<T>(c_phi_evo(k, l, n))
						* point_to_ellipse_series::series_pow<T>(abs_sin_psi, l)
						* point_to_ellipse_series::series_pow<T>(rho_ae2, k)
						* point_to_ellipse_series::series_pow<T>(b_a, n);
	return d;
}

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
inline T phi_evo_dense(int K,
					   const T& abs_sin_psi, const T& rho_ae2, const T& b_a) {
	T d(0);
	for (int k = 1; k <= K; ++k) {
		T cm(0);
		const int s = (k + 1) % 2;
		const int r = (k - 1) / 2;
		for (int l = 0; l <= r; ++l) {
			for (int n = 0; n <= r; ++n)
				cm = cm + point_to_ellipse_series::series_coeff<T>(d_phi_evo2(k, l, n))
						  * point_to_ellipse_series::series_pow<T>(abs_sin_psi, 2 * l)
						  * point_to_ellipse_series::series_pow<T>(b_a, 2 * n);
		}
		d = d + cm * point_to_ellipse_series::series_pow<T>(rho_ae2, k)
				* point_to_ellipse_series::series_pow<T>(abs_sin_psi, s)
				* point_to_ellipse_series::series_pow<T>(b_a, s + 1);
	}
	return d;
}

template<typename T>
inline T sin_phi_evo_sparse(int L, int K,
							const T& abs_sin_psi,
							const T& rho_ae2, const T& b_a) {
	T d(0);
	for (int l = 0; l <= L; ++l)
		for (int k = l; k <= K; ++k)
			for (int n = 0; n <= k; ++n)
				d = d + point_to_ellipse_series::series_coeff<T>(c_sin_phi_evo(k, l, n))
						* point_to_ellipse_series::series_pow<T>(abs_sin_psi, l)
						* point_to_ellipse_series::series_pow<T>(rho_ae2, k)
						* point_to_ellipse_series::series_pow<T>(b_a, n);
	return d;
}


/** Inside-evolute series for sin(phi) in sin-powers (dense).
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param K     rho_ae2 power series truncation order.
 * @param abs_sin_psi   Value/expression for |sin(psi)|.
 * @param rho_ae2   Value/expression for rho/(a*e²).
 * @param b_a       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T sin_phi_evo_dense(int K,
						   const T& abs_sin_psi, const T& rho_ae2, const T& b_a) {
	T s(0);
	for (int k = 2; k <= K; ++k) {
		T cm(0);
		const int prtm = k % 2;
		for (int l = 0; l <= k / 2; ++l) {
			for (int n = 1; n <= k / 2; ++n) {
				cm = cm
					 + point_to_ellipse_series::series_coeff<T>(d_sin_phi_evo2(k, l, n))
					   * point_to_ellipse_series::series_pow<T>(b_a, 2 * n)
					   * point_to_ellipse_series::series_pow<T>(abs_sin_psi, 2 * l);
			}
		}
		s = s + cm
				* point_to_ellipse_series::series_pow<T>(b_a, prtm)
				* point_to_ellipse_series::series_pow<T>(abs_sin_psi, prtm)
				* point_to_ellipse_series::series_pow<T>(rho_ae2, k);
	}
	return s;
}

template<typename T>
inline T cos_phi_evo_sparse(int L, int K,
							const T& abs_sin_psi, const T& rho_ae2, const T& b_a) {
	T d(0);
	for (int l = 0; l <= L; ++l)
		for (int k = l; k <= K; ++k)
			for (int n = 1; n <= k; ++n)
				d = d + point_to_ellipse_series::series_coeff<T>(c_cos_phi_evo(k, l, n))
						* point_to_ellipse_series::series_pow<T>(abs_sin_psi, l)
						* point_to_ellipse_series::series_pow<T>(rho_ae2, k)
						* point_to_ellipse_series::series_pow<T>(b_a, n);
	return d;
}


/** Inside-evolute series for cos(phi) / |cos(psi)| in sin-powers (dense).
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho_ae2 power series truncation order.
 * @param abs_sin_psi   Value/expression for |sin(psi)|.
 * @param rho_ae2   Value/expression for rho/(a*e²).
 * @param b_a       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T cos_phi_evo_dense(int K,
						   const T& abs_sin_psi,
						   const T& rho_ae2,
						   const T& b_a) {
	T s(0);
	for (int k = 1; k <= K; ++k) {
		T cm(0);

		const int p = (k - 1) % 2;
		const int q  = (k - 1) / 2;

		for (int l = 0; l <= q; ++l) {
			for (int n = p; n <= k / 2; ++n) {
				cm = cm
					 + point_to_ellipse_series::series_coeff<T>(d_cos_phi_evo_m(k, l, n))
					   * point_to_ellipse_series::series_pow<T>(b_a, 2 * n)
					   * point_to_ellipse_series::series_pow<T>(abs_sin_psi, 2 * l);
			}
		}
		s = s + cm * point_to_ellipse_series::series_pow<T>(b_a, 1 - p)
				* point_to_ellipse_series::series_pow<T>(abs_sin_psi, p)
				* point_to_ellipse_series::series_pow<T>(rho_ae2, k);
	}

	return s;
}

/** Inside-evolute series for h/a - rho/a*sin(psi) in sin-powers.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param L     sin power series truncation order.
 * @param K     rho_ae2 power series truncation order.
 * @param abs_sin_psi   Value/expression for |sin(psi)|.
 * @param rho_ae2   Value/expression for rho/(a*e²).
 * @param b_a       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T h_evo_sparse(int L, int K,
					  const T& abs_sin_psi, const T& rho_ae2, const T& b_a) {
	T s(0);
	for (int l = 0; l <= L; ++l)
		for (int k = l; k <= K; ++k)
			for (int n = 1; n <= k + 1; ++n)
				s = s + point_to_ellipse_series::series_coeff<T>(c_h_evo2(k, l, n))
						* point_to_ellipse_series::series_pow<T>(b_a, n)
						* point_to_ellipse_series::series_pow<T>(rho_ae2, k)
						* point_to_ellipse_series::series_pow<T>(abs_sin_psi, l);
	return s;
}

/** Inside-evolute series for (h - b + rho * |sin(psi)|) / a in sin-powers (dense).
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
inline T h_evo_dense(int K,
					 const T& sin_psi_v, const T& rho_ae2_v, const T& b_a_v) {
	T s(0);

	for (int k = 2; k <= K; ++k) {
		T cm(0);
		const int prtm = k % 2;

		for (int l = 0; l <= k / 2; ++l) {
			for (int n = 0; n <= k / 2; ++n) {
				cm = cm + point_to_ellipse_series::series_coeff<T>(d_h_evo3(k, l, n))
						  * point_to_ellipse_series::series_pow<T>(b_a_v, 2 * n)
						  * point_to_ellipse_series::series_pow<T>(sin_psi_v, 2 * l);
			}
		}
		s = s + cm * point_to_ellipse_series::series_pow<T>(b_a_v, 1+prtm)
				* point_to_ellipse_series::series_pow<T>(sin_psi_v, prtm)
				* point_to_ellipse_series::series_pow<T>(rho_ae2_v, k);
	}
	return s;
}

