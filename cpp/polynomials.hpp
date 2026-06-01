#pragma once

#include <cassert>
#include <functional>
#include <algorithm>
#include <map>
#include <tuple>

#include "cache.hpp"
#include "util.hpp"
#include "lexpr.hpp"

namespace point_to_ellipse_series {

/**
 * Fourier multiple-angle cosine series coefficient from sin-power series.
 * GMP-only — returns mpq_class directly.
 *
 * @param n      Sin-multiple.
 * @param k      rho power.
 * @param l      e² power.
 * @param d_nkl  Sin-power series coefficient function: (n, k, l) → mpq_class.
 * @param n_min  Lowest sin-power.
 * @param k_pp   Maximum e² power offset in sin-power series.
 * @return       Fourier multiple-angle cos series coefficient.
 */
inline mpq_class
sin_pow_to_cos_mul(int n, int k, int l, int n_min, int k_pp,
				   const std::function<mpq_class(int, int, int)> &d_nkl) {
	assert(n >= 0 && k >= 0 && l >= 0);
	mpq_class c_nkl(0);
	for (int i = std::max({n, n_min, l - k - k_pp}); i <= l; ++i) {
		mpz_class binom;
		mpz_bin_uiui(binom.get_mpz_t(), (unsigned long)(2 * i), (unsigned long)(i - n));
		mpz_class denom;
		mpz_ui_pow_ui(denom.get_mpz_t(), 2, (unsigned long)(2 * i));
		mpq_class term = d_nkl(i, k, l) * mpq_class(binom, denom);
		term.canonicalize();
		c_nkl += term;
	}
	c_nkl *= mpq_class((1 + !!n) * powm1(n));
	c_nkl.canonicalize();
	return c_nkl;
}

// ---------------------------------------------------------------------------
// GMP-only polynomial type for a_j variables (index → exponent).
// ---------------------------------------------------------------------------

// Monomial: maps variable index j to its exponent.
using AMonomial = std::map<int, int>;
// Polynomial: maps monomial to its rational coefficient.
using APoly = std::map<AMonomial, mpq_class>;

// Add (scale * a_j * q) into p in-place.
inline static void apoly_add_scaled_mul_var(APoly &p, const APoly &q, int j,
											const mpq_class &scale) {
	for (const auto &[mono, coeff] : q) {
		AMonomial new_mono = mono;
		new_mono[j]++;
		mpq_class c = scale * coeff;
		c.canonicalize();
		p[new_mono] += c;
	}
}

// Divide every monomial by a_0 (asserts a_0 exponent >= 1 in each term).
inline static APoly apoly_div_a0(const APoly &p) {
	APoly result;
	for (const auto &[mono, coeff] : p) {
		if (coeff == mpq_class(0)) continue;
		AMonomial new_mono = mono;
		auto it = new_mono.find(0);
		assert(it != new_mono.end() && it->second >= 1);
		if (--it->second == 0)
			new_mono.erase(it);
		result[new_mono] += coeff;
	}
	return result;
}

// Compute the n-th coefficient polynomial of (sum_j a_j x^j)^i using GMP arithmetic.
// Returns a polynomial in indexed a_j variables with mpq_class coefficients.
inline APoly ordinary_potential_polynomial2(int n, int i) {
	assert(n >= 0 && i >= 0);
	static UintsCache<APoly> cache;
	if (auto *ret = cache.get((uint) n, (uint) i))
		return *ret;

	APoly A;
	if (n == 0) {
		AMonomial mono;
		if (i > 0) mono[0] = i;
		A[mono] = mpq_class(1);
	} else {
		APoly clj;
		for (int k = 1; k <= n; ++k) {
			long factor = (long) k * i - n + k;
			if (factor == 0) continue;
			apoly_add_scaled_mul_var(clj, ordinary_potential_polynomial2(n - k, i), k,
									 mpq_class(factor));
		}
		APoly div = apoly_div_a0(clj);
		mpq_class inv_n(1, n);
		for (auto &[mono, coeff] : div) {
			coeff *= inv_n;
			coeff.canonicalize();
		}
		A = std::move(div);
	}

	return cache.insert(A, (uint) n, (uint) i);
}

// ---------------------------------------------------------------------------
// Dense polynomial in e² (index = power of e²).
// ---------------------------------------------------------------------------

using E2Poly = std::vector<mpq_class>;

inline E2Poly e2poly_add(E2Poly a, const E2Poly& b) {
	if (b.size() > a.size()) a.resize(b.size(), mpq_class(0));
	for (size_t i = 0; i < b.size(); ++i) {
		a[i] += b[i];
		a[i].canonicalize();
	}
	return a;
}

inline E2Poly e2poly_mul(const E2Poly& a, const E2Poly& b) {
	if (a.empty() || b.empty()) return {};
	E2Poly result(a.size() + b.size() - 1, mpq_class(0));
	for (size_t i = 0; i < a.size(); ++i) {
		if (a[i] == 0) continue;
		for (size_t j = 0; j < b.size(); ++j) {
			if (b[j] == 0) continue;
			result[i + j] += a[i] * b[j];
			result[i + j].canonicalize();
		}
	}
	return result;
}

// Partial ordinary Bell polynomial B_{k,i}(a_{j,1}, a_{j,2}, ...) as LExpr.
// Variables: lexpr_var(j, m) for inner index m = 1, 2, ..., k-i+1.
inline LExpr partial_bell_lexpr(int k, int i, int j) {
	static std::map<std::tuple<int,int,int>, LExpr> cache;
	auto key = std::make_tuple(k, i, j);
	auto it = cache.find(key);
	if (it != cache.end()) return it->second;

	LExpr B;
	if (i == 0) {
		B = (k == 0) ? lexpr_const(mpq_class(1)) : LExpr{};
	} else {
		for (int m = 1; m <= k - i + 1; ++m)
			B += lexpr_var(j, m) * partial_bell_lexpr(k - m, i - 1, j);
	}
	return cache[key] = B;
}

} // namespace point_to_ellipse_series
