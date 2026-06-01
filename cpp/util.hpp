#pragma once

#include <unordered_map>
#include <gmpxx.h>

#include "coefficients.hpp"

namespace point_to_ellipse_series {

inline std::string mpz_to_str(const mpz_class& x) { return x.get_str(10); }

/** Compute the n:th rising factorial of (k/2).
 *
 * @param k
 * @param n
 * @return
 */
inline mpq_class rf_half(unsigned long k, unsigned long n) {
	mpz_class r(1), s(k);
	for (unsigned long i = 0; i < n; ++i) { r *= s; s += 2; }
	mpz_class denom;
	mpz_ui_pow_ui(denom.get_mpz_t(), 2, n);
	mpq_class result(r, denom);
	result.canonicalize();
	return result;
}

/** (-1)^n
 *
 * See
 * https://stackoverflow.com/questions/29110752/what-is-the-correct-way-to-obtain-1n
 *
 * @param n
 * @return
 */
inline long powm1(long n) {
	return 1L - ((n & 1L) << 1);
}

/** Generalised integer binomial C(n, k) for arbitrary integer n and any k.
 *
 * Returns 0 for k < 0. For n < 0 uses C(n,k) = (-1)^k * C(k-n-1, k).
 *
 * @param n Integer (may be negative).
 * @param k Integer.
 * @return
 */
inline mpq_class int_bin(long n, long k) {
	if (k < 0) return {0};
	if (n >= 0) {
		if (k > n) return {0};
		mpz_class result;
		mpz_bin_uiui(result.get_mpz_t(), (unsigned long) n, (unsigned long) k);
		return mpq_class(result);
	}
	// n < 0: C(n,k) = (-1)^k * C(k-n-1, k)
	mpz_class result;
	mpz_bin_uiui(result.get_mpz_t(), (unsigned long) (k - n - 1), (unsigned long) k);
	return mpq_class(mpz_class(powm1(k)) * result);
}

/** Compute the binomial coefficient for rational n and integer k.
 *
 * Computes (n-k+1)*(n-k+2)*...*n / k!  Returns 0 for k < 0.
 *
 * @param n Rational number.
 * @param k Integer.
 * @return
 */
inline mpq_class binomial_rational(const mpq_class& n, long k) {
	if (k < 0) return {0};
	if (k == 0) return {1};
	mpq_class result(1);
	mpq_class d = n - mpq_class(k - 1);
	for (long i = 0; i < k; ++i) {
		result *= d;
		d += 1;
	}
	mpz_class fact;
	mpz_fac_ui(fact.get_mpz_t(), (unsigned long) k);
	result /= mpq_class(fact);
	result.canonicalize();
	return result;
}

/** Euler secant number E_{2n}.
 *
 * Mirrors Python E2(n). Uses a static cache for efficiency.
 * Defined by: E2(0)=1, E2(n) = sum_{j=0}^{n-1} (-1)^(n-j+1) * C(2n,2j) * E2(j)
 *
 * @param n Non-negative integer.
 * @return E_{2n} as mpz_class.
 */
inline mpz_class E2(int n) {
	static std::unordered_map<int, mpz_class> cache;
	auto it = cache.find(n);
	if (it != cache.end())
		return it->second;

	if (n == 0)
		return mpz_class(1);

	mpz_class result(0);
	for (int j = 0; j < n; ++j) {
		mpz_class binom;
		mpz_bin_uiui(binom.get_mpz_t(), (unsigned long)(2 * n), (unsigned long)(2 * j));
		result += mpz_class(powm1(n - j + 1)) * binom * E2(j);
	}

	cache[n] = result;
	return result;
}

} // namespace point_to_ellipse_series
