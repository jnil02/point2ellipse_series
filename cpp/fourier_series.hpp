#pragma once

/*
 * Symbolic series expansions for the point-to-ellipse relation.
 */

#include "coefficients.hpp"
#include "polynomials.hpp"
#include "series_traits.hpp"

using point_to_ellipse_series::d_phi;
using point_to_ellipse_series::d_phi2;
using point_to_ellipse_series::c_phi;
using point_to_ellipse_series::d_sin;
using point_to_ellipse_series::c_sin;
using point_to_ellipse_series::d_cos;
using point_to_ellipse_series::c_cos;
using point_to_ellipse_series::d_h;
using point_to_ellipse_series::c_h;

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
				d = d + point_to_ellipse_series::series_coeff<T>(d_phi(n, k, l))
						* point_to_ellipse_series::series_pow<T>(e2_v, l)
						* point_to_ellipse_series::series_pow<T>(varrho_v, k)
						* point_to_ellipse_series::series_pow<T>(sin_psi_v, 2 * n);
	return d;
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
				d = d + point_to_ellipse_series::series_coeff<T>(d_phi2(n, k, l))
						* point_to_ellipse_series::series_pow<T>(e2_v, l)
						* point_to_ellipse_series::series_pow<T>(varrho_v, k)
						* point_to_ellipse_series::series_pow<T>(sin_psi_v, 2 * n + 1);
	return d;
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
				d = d + point_to_ellipse_series::series_coeff<T>(c_phi(n, k, l))
						* point_to_ellipse_series::series_pow<T>(e2_v, l)
						* point_to_ellipse_series::series_pow<T>(varrho_v, k)
						* point_to_ellipse_series::series_sin_mul<T>(psi_v, n);
	return d;
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
				d = d + point_to_ellipse_series::series_coeff<T>(d_sin(n, k, l))
						* point_to_ellipse_series::series_pow<T>(e2_v, l)
						* point_to_ellipse_series::series_pow<T>(varrho_v, k)
						* point_to_ellipse_series::series_pow<T>(sin_psi_v, 2 * n);
	return d;
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
				d = d + point_to_ellipse_series::series_coeff<T>(c_sin(n, k, l))
						* point_to_ellipse_series::series_pow<T>(e2_v, l)
						* point_to_ellipse_series::series_pow<T>(varrho_v, k)
						* point_to_ellipse_series::series_cos_mul<T>(psi_v, n);
	return d;
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
				d = d + point_to_ellipse_series::series_coeff<T>(d_cos(n, k, l))
						* point_to_ellipse_series::series_pow<T>(e2_v, l)
						* point_to_ellipse_series::series_pow<T>(varrho_v, k)
						* point_to_ellipse_series::series_pow<T>(sin_psi_v, 2 * n);
	return d;
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
				d = d + point_to_ellipse_series::series_coeff<T>(c_cos(n, k, l))
						* point_to_ellipse_series::series_pow<T>(e2_v, l)
						* point_to_ellipse_series::series_pow<T>(varrho_v, k)
						* point_to_ellipse_series::series_cos_mul<T>(psi_v, n);
	return d;
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
				d = d + point_to_ellipse_series::series_coeff<T>(d_h(n, k, l))
						* point_to_ellipse_series::series_pow<T>(e2_v, l)
						* point_to_ellipse_series::series_pow<T>(varrho_v, k)
						* point_to_ellipse_series::series_pow<T>(sin_psi_v, 2 * n);
	return d;
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
				d = d + point_to_ellipse_series::series_coeff<T>(c_h(n, k, l))
						* point_to_ellipse_series::series_pow<T>(e2_v, l)
						* point_to_ellipse_series::series_pow<T>(varrho_v, k)
						* point_to_ellipse_series::series_cos_mul<T>(psi_v, n);
	return d;
}
