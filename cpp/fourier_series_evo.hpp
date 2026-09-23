#pragma once

/*
 * Symbolic inside-evolute series expansions for the point-to-ellipse relation.
 */

#include "fourier_series.hpp"
#include "coefficients_evo.hpp"

namespace point_to_ellipse_series {

/** Inside-evolute series for (phi - sgn*pi/2) / (sgn*|cos(psi)|) in sin-powers (sparse).
 *
 * Matches Python fourier_series.phi_evo_sparse.
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
inline T phi_evo_sparse(int L, int K,
						const T &abs_sin_psi,
						const T &rho_ae2,
						const T &b_a) {
	T s(0);
	for (int l = 0; l <= L; ++l)
		for (int k = l + 1; k <= K; ++k)
			for (int n = 1; n <= k; ++n)
				s += to<T>(c_phi_evo(k, l, n))
					 * ipow<T>(b_a, n)
					 * ipow<T>(abs_sin_psi, l)
					 * ipow<T>(rho_ae2, k);
	return s;
}

/** Inside-evolute series for (phi - sgn*pi/2) / (sgn*|cos(psi)|) in sin-powers (dense).
 *
 * Matches Python fourier_series.phi_evo_dense.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param K     rho_ae2 power series truncation order.
 * @param abs_sin_psi   Value/expression for |sin(psi)|.
 * @param rho_ae2   Value/expression for rho/(a*e²).
 * @param b_a       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T phi_evo_dense(int K,
					   const T &abs_sin_psi,
					   const T &rho_ae2,
					   const T &b_a) {
	T s(0);
	for (int k = 1; k <= K; ++k) {

		T cm(0);
		const int p = (k + 1) % 2;
		const int q = (k - 1) / 2;

		for (int l = 0; l <= q; ++l) {
			for (int n = 0; n <= q; ++n)
				cm += to<T>(d_phi_evo(k, l, n))
					  * ipow<T>(b_a, 2 * n)
					  * ipow<T>(abs_sin_psi, 2 * l);
		}
		s += cm * ipow<T>(b_a, p + 1)
			 * ipow<T>(abs_sin_psi, p)
			 * ipow<T>(rho_ae2, k);
	}
	return s;
}

/** Inside-evolute series for sin(phi) in sin-powers (sparse).
 *
 * Matches Python fourier_series.sin_phi_evo_sparse.
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
inline T sin_phi_evo_sparse(int L, int K,
							const T &abs_sin_psi,
							const T &rho_ae2,
							const T &b_a) {
	T s(0);
	for (int l = 0; l <= L; ++l)
		for (int k = l; k <= K; ++k)
			for (int n = 0; n <= k; ++n)
				s += to<T>(c_sin_phi_evo(k, l, n))
					 * ipow<T>(b_a, n)
					 * ipow<T>(abs_sin_psi, l)
					 * ipow<T>(rho_ae2, k);
	return s;
}


/** Inside-evolute series for sin(phi) in sin-powers (dense).
 *
 * Matches Python fourier_series.sin_phi_evo_dense.
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
						   const T &abs_sin_psi,
						   const T &rho_ae2,
						   const T &b_a) {
	T s(0);
	for (int k = 2; k <= K; ++k) {
		T cm(0);
		const int p = k % 2;
		const int q = k / 2;
		for (int l = 0; l <= q; ++l) {
			for (int n = 1; n <= q; ++n) {
				cm += to<T>(d_sin_phi_evo(k, l, n))
					  * ipow<T>(b_a, 2 * n)
					  * ipow<T>(abs_sin_psi, 2 * l);
			}
		}
		s += cm * ipow<T>(b_a, p)
			 * ipow<T>(abs_sin_psi, p)
			 * ipow<T>(rho_ae2, k);
	}
	return s;
}

/** Inside-evolute series for cos(phi) / |cos(psi)| in sin-powers (sparse).
 *
 * Matches Python fourier_series.cos_phi_evo_sparse.
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
inline T cos_phi_evo_sparse(int L, int K,
							const T &abs_sin_psi,
							const T &rho_ae2,
							const T &b_a) {
	T d(0);
	for (int l = 0; l <= L; ++l)
		for (int k = l; k <= K; ++k)
			for (int n = 1; n <= k; ++n)
				d += to<T>(c_cos_phi_evo(k, l, n))
					 * ipow<T>(b_a, n)
					 * ipow<T>(abs_sin_psi, l)
					 * ipow<T>(rho_ae2, k);
	return d;
}


/** Inside-evolute series for cos(phi) / |cos(psi)| in sin-powers (dense).
 *
 * Matches Python fourier_series.cos_phi_evo_dense.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param K     rho_ae2 power series truncation order.
 * @param abs_sin_psi   Value/expression for |sin(psi)|.
 * @param rho_ae2   Value/expression for rho/(a*e²).
 * @param b_a       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T cos_phi_evo_dense(int K,
						   const T &abs_sin_psi,
						   const T &rho_ae2,
						   const T &b_a) {
	T s(0);
	for (int k = 1; k <= K; ++k) {
		T cm(0);
		const int p = (k - 1) % 2;
		const int q = (k - 1) / 2;
		for (int l = 0; l <= q; ++l) {
			for (int n = p; n <= k / 2; ++n) {
				cm += to<T>(d_cos_phi_evo(k, l, n))
					  * ipow<T>(b_a, 2 * n)
					  * ipow<T>(abs_sin_psi, 2 * l);
			}
		}
		s += cm * ipow<T>(b_a, 1 - p)
			 * ipow<T>(abs_sin_psi, p)
			 * ipow<T>(rho_ae2, k);
	}

	return s;
}

/** Inside-evolute series for h/a - rho/a*sin(psi) in sin-powers.
 *
 * Matches Python fourier_series.h_evo_sparse.
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
					  const T &abs_sin_psi,
					  const T &rho_ae2,
					  const T &b_a) {
	T s(0);
	for (int l = 0; l <= L; ++l)
		for (int k = l; k <= K; ++k)
			for (int n = 1; n <= k + 1; ++n)
				s += to<T>(c_h_evo(k, l, n))
					 * ipow<T>(b_a, n)
					 * ipow<T>(abs_sin_psi, l)
					 * ipow<T>(rho_ae2, k);
	return s;
}

/** Inside-evolute series for (h + b - rho * |sin(psi)|) / a in sin-powers (dense).
 *
 * Matches Python fourier_series.h_evo_dense.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param K     rho_ae2 power series truncation order.
 * @param abs_sin_psi   Value/expression for |sin(psi)|.
 * @param rho_ae2   Value/expression for rho/(a*e²).
 * @param b_a       Value/expression for b/a.
 * @return Series result as type T.
 */
template<typename T>
inline T h_evo_dense(int K,
					 const T &abs_sin_psi,
					 const T &rho_ae2,
					 const T &b_a) {
	T s(0);

	for (int k = 2; k <= K; ++k) {
		T cm(0);
		const int p = k % 2;
		const int q = k / 2;
		for (int l = 0; l <= q; ++l) {
			for (int n = 0; n <= q; ++n) {
				cm += to<T>(d_h_evo(k, l, n))
					  * ipow<T>(b_a, 2 * n)
					  * ipow<T>(abs_sin_psi, 2 * l);
			}
		}
		s += cm * ipow<T>(b_a, 1 + p)
			 * ipow<T>(abs_sin_psi, p)
			 * ipow<T>(rho_ae2, k);
	}
	return s;
}

} // namespace point_to_ellipse_series
