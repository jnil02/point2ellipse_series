#pragma once

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/pow.h>
#include <symengine/ntheory.h>
#include <symengine/number.h>
#include <functional>
#include <algorithm>

#include "cache.hpp"
#include "util.hpp"

using SymEngine::Expression;
using SymEngine::integer;
using SymEngine::binomial;
using SymEngine::Integer;
using SymEngine::rational;
using SymEngine::symbol;
using SymEngine::integer;
using SymEngine::div;
using SymEngine::factorial;

/**
 * Fourier multiple-angle cosine series coefficient from sin-power series.
 *
 * @param n      Sin-multiple.
 * @param k      rho power.
 * @param l      e² power.
 * @param d_nkl  Sin-power series coefficient function: (n, k, l) → Expression.
 * @param n_min  Lowest sin-power.
 * @param k_pp   Maximum e² power offset in sin-power series.
 * @return       Fourier multiple-angle cos series coefficient.
 */
inline Expression
sin_pow_to_cos_mul(int n, int k, int l, int n_min, int k_pp,
				   const std::function<Expression(int, int, int)> &d_nkl) {
	assert(n >= 0 && k >= 0 && l >= 0);

	Expression c_nkl(0);
	for (int i = std::max({n, n_min, l - k - k_pp}); i <= l; ++i)
		c_nkl = c_nkl + d_nkl(i, k, l) * rational(1, 1L << (2 * i))
				* binomial(Integer(2 * i), (unsigned long)(i - n));
	// Common factor 2*(-1)^n / (n == 0 ? 2 : 1)
	c_nkl = c_nkl * ((1 + !!n) * powm1(n));

	return c_nkl;
}

/**
 * Compute the partial ordinary Bell polynomial recursively.
 *
 * @param k Power index of the polynomial.
 * @param i Sum index of the polynomial.
 * @param a Base name for the variables (e.g., "a" → a_1, a_2, ...).
 * @return SymEngine Expression for the partial ordinary Bell polynomial.
 */
inline Expression partial_ordinary_bell_polynomial(int k, int i, const std::string& a) {

	// Look in cache for result.
	static std::map<std::tuple<int,std::string>, Expression> cache;
	auto key = std::make_tuple(cantor_pairing_two(k, i), a);
	auto it = cache.find(key);
	if (it != cache.end())
		return it->second;

	// Recursive computation of ordinary Bell polynomial.
	Expression B;
	if (i == 0) {
		B = {(k == 0) ? integer(1) : integer(0)};
	} else {
		B = {integer(0)};
		for (int j = 1; j <= k - i + 1; ++j)
			B += Expression(symbol(a + "_" + std::to_string(j)))
				 * partial_ordinary_bell_polynomial(k - j, i - 1, a);
	}

	// Expand and store in cache before returning.
	cache[key] = expand(B);
	return B;
}

/**
 * Compute the n-th coefficient polynomial of the i-th power of an infinite power series.
 *
 * Given:
 *   p = sum_{n=0}^\infty a_n x^n  (with a_0 != 0)
 *   p^i = sum_{n=0}^\infty b_{n,i} x^n
 *
 * This computes the polynomial P_{n,i}(a_0, ..., a_n) for b_{n,i}.
 *
 * @param n The index of the coefficient.
 * @param i The power of the original series.
 * @param a The base symbol name of the coefficients (e.g., "a" → a_0, a_1, ...).
 * @return SymEngine Expression for b_{n,i}.
 */
inline Expression ordinary_potential_polynomial(int n, int i, const std::string& a) {

	// Look in cache for result.
	static std::map<std::tuple<int,std::string>, Expression> cache;
	auto key = std::make_tuple(cantor_pairing_two(n, i), a);
	auto it = cache.find(key);
	if (it != cache.end())
		return it->second;

	// Define x0 = a_0.
	Expression x0(symbol(a + "_0"));

	// Recursive computation of ordinary potential polynomial.
	Expression A;
	if (n == 0) {
		A = expand(pow(x0, i));
	} else {
		Expression clj(integer(0));
		for (int k = 1; k <= n; ++k)
			clj = clj + integer(k * i - n + k)
					* Expression(a + "_" + std::to_string(k)) *
					ordinary_potential_polynomial(n - k, i, a);

		A = clj / (integer(n) * x0);
	}

	// Expand and store in cache before returning.
	cache[key] = expand(A);
	return A;
}

// "Sigma" polynomial giving one common component of the sin(phi)/sin(psi) and the cos(phi)/cos(psi) expansions.
inline Expression sigma(int J, const Expression &delta) {
	Expression d(0);
	for (int j = 0; j < J; j += 2) {
		int j2 = j / 2;
		d += Expression(div(integer(powm1(j2)), factorial(j))) * pow(delta, j2);
	}
	return d;
}

// "Tau" polynomial giving one common component of the sin(phi)/sin(psi) and the cos(phi)/cos(psi) expansions.
inline Expression
tau(int J, const Expression &omega, const Expression &delta) {
	Expression d(0);
	for (int j = 1; j < J; j += 2) {
		int j2 = (j - 1) / 2;
		d += Expression(div(integer(powm1(j2)), factorial(j))) * pow(delta, j2);
	}
	return omega * d;
}
