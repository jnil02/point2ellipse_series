#pragma once

/*
 * Symbolic series expansions for the point-to-ellipse relation.
 */

#include "coefficients.hpp"
#include "series_traits.hpp"

#include "detail/polynomials.hpp"

namespace point_to_ellipse_series {

/** Series expansion of (phi - psi) / (sin(psi) * cos(psi)) in sin-powers.
 *
 * Matches Python expansions.phi_in_sin_pow.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for a/rho.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T phi_in_sin_pow(int N, int K,
						const T &sin_psi_v, const T &varrho_v, const T &e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(k, n + 1); l <= k + n; ++l)
				d += to<T>(d_phi(n, k, l))
					 * ipow<T>(e2_v, l)
					 * ipow<T>(varrho_v, k)
					 * ipow<T>(sin_psi_v, 2 * n);
	return d;
}

/** Series expansion of (phi - psi) in sin-powers.
 *
 * Matches Python expansions.phi_in_sin_pow2.
 *
 * Note, this series has poor convergence and is only implemented to demonstrate
 * this. It should not be used in practice.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for a/rho.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T phi_in_sin_pow2(int N, int K,
						 const T &sin_psi_v, const T &varrho_v, const T &e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = k; l <= k + n; ++l)
				d += to<T>(d_phi2(n, k, l))
					 * ipow<T>(e2_v, l)
					 * ipow<T>(varrho_v, k)
					 * ipow<T>(sin_psi_v, 2 * n + 1);
	return d;
}

/** Series expansion of (phi - psi) in sin multiples.
 *
 * Matches Python expansions.phi_in_sin_mul.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin multiple series truncation order.
 * @param K     rho power series truncation order.
 * @param L     e2 power series truncation order.
 * @param psi_v       Value/expression for psi.
 * @param varrho_v    Value/expression for a/rho.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T phi_in_sin_mul(int N, int K, int L,
						const T &psi_v, const T &varrho_v, const T &e2_v) {
	T d(0);
	for (int n = 1; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(n, k); l <= L; ++l)
				d += to<T>(c_phi(n, k, l))
					 * ipow<T>(e2_v, l)
					 * ipow<T>(varrho_v, k)
					 * isin_mul<T>(psi_v, n);
	return d;
}

/** Series expansion of sin(phi) / sin(psi) - 1 in sin-powers.
 *
 * Matches Python expansions.sin_phi_in_sin_pow.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for a/rho.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T sin_phi_in_sin_pow(int N, int K,
							const T &sin_psi_v, const T &varrho_v,
							const T &e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(k, n); l <= n + k; ++l)
				d += to<T>(d_sin(n, k, l))
					 * ipow<T>(e2_v, l)
					 * ipow<T>(varrho_v, k)
					 * ipow<T>(sin_psi_v, 2 * n);
	return d;
}

/** Series expansion of sin(phi) / sin(psi) - 1 in cos multiples.
 *
 * Matches Python expansions.sin_phi_in_cos_mul.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin multiple series truncation order.
 * @param K     rho power series truncation order.
 * @param L     e2 power series truncation order.
 * @param psi_v       Value/expression for psi.
 * @param varrho_v    Value/expression for a/rho.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T sin_phi_in_cos_mul(int N, int K, int L,
							const T &psi_v, const T &varrho_v, const T &e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(n, k); l <= L; ++l)
				d += to<T>(c_sin(n, k, l))
					 * ipow<T>(e2_v, l)
					 * ipow<T>(varrho_v, k)
					 * icos_mul<T>(psi_v, n);
	return d;
}

/** Series expansion of cos(phi) / cos(psi) - 1 in sin-powers.
 *
 * Matches Python expansions.cos_phi_in_sin_pow.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for a/rho.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T cos_phi_in_sin_pow(int N, int K,
							const T &sin_psi_v, const T &varrho_v,
							const T &e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(k, n); l < n + k; ++l)
				d += to<T>(d_cos(n, k, l))
					 * ipow<T>(e2_v, l)
					 * ipow<T>(varrho_v, k)
					 * ipow<T>(sin_psi_v, 2 * n);
	return d;
}

/** Series expansion of cos(phi) / cos(psi) - 1 in cos multiples.
 *
 * Matches Python expansions.cos_phi_in_cos_mul.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin multiple series truncation order.
 * @param K     rho power series truncation order.
 * @param L     e2 power series truncation order.
 * @param psi_v       Value/expression for psi.
 * @param varrho_v    Value/expression for a/rho.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T cos_phi_in_cos_mul(int N, int K, int L,
							const T &psi_v, const T &varrho_v, const T &e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 1; k <= K; ++k)
			for (int l = std::max(n, k); l <= L; ++l)
				d += to<T>(c_cos(n, k, l))
					 * ipow<T>(e2_v, l)
					 * ipow<T>(varrho_v, k)
					 * icos_mul<T>(psi_v, n);
	return d;
}

/** Series expansion of (h + a - rho) / a in sin-powers.
 *
 * Matches Python expansions.h_in_sin_pow.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin power series truncation order.
 * @param K     rho power series truncation order.
 * @param sin_psi_v   Value/expression for sin(psi).
 * @param varrho_v    Value/expression for a/rho.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T h_in_sin_pow(int N, int K,
					  const T &sin_psi_v, const T &varrho_v, const T &e2_v) {
	T d(0);
	for (int n = 1; n <= N; ++n)
		for (int k = 0; k <= K; ++k)
			for (int l = std::max(k + 1, n); l <= n + k; ++l)
				d += to<T>(d_h(n, k, l))
					 * ipow<T>(e2_v, l)
					 * ipow<T>(varrho_v, k)
					 * ipow<T>(sin_psi_v, 2 * n);
	return d;
}

/** Series expansion of (h + a - rho) / a in cos multiples.
 *
 * Matches Python expansions.h_in_cos_mul.
 *
 * @tparam T    Value type: SymEngine::Expression for symbolic, mpfr::mpreal for numeric.
 * @param N     sin multiples series truncation order.
 * @param K     rho power series truncation order.
 * @param L     e2 power series truncation order.
 * @param psi_v       Value/expression for psi.
 * @param varrho_v    Value/expression for a/rho.
 * @param e2_v        Value/expression for e².
 * @return Series result as type T.
 */
template<typename T>
inline T h_in_cos_mul(int N, int K, int L,
					  const T &psi_v, const T &varrho_v, const T &e2_v) {
	T d(0);
	for (int n = 0; n <= N; ++n)
		for (int k = 0; k <= K; ++k)
			for (int l = std::max(n, k + 1); l <= L; ++l)
				d += to<T>(c_h(n, k, l))
					 * ipow<T>(e2_v, l)
					 * ipow<T>(varrho_v, k)
					 * icos_mul<T>(psi_v, n);
	return d;
}

} // namespace point_to_ellipse_series