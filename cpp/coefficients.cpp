
#include <symengine/expression.h>
#include <symengine/ntheory.h>
#include <symengine/rational.h>
#include <symengine/symbol.h>

#include "series.hpp"
#include "series_substitution.hpp"
#include "cache.hpp"
#include "util.hpp"

#include "coefficients.hpp"

namespace point_to_ellipse_series {

using SymEngine::Expression;
using SymEngine::Rational;
using SymEngine::rational;
using SymEngine::binomial;
using SymEngine::rcp_static_cast;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::Integer;
using SymEngine::integer;
using SymEngine::factorial;
using SymEngine::add;
using SymEngine::rational;

using uint = unsigned int;

rc d_phi(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d = Expression(0);
	for (int r = 0; r <= k - 1; ++r)
		for (int m = 0; m <= k - 1 - r; ++m)
			for (int q = 0; q <= r / 2; ++q)
				for (int p = 0; p <= m / 2; ++p)
					for (int t = l - n - k + m + r - p - q; t <= std::min(
							{l - k, r - q, l - n - k + m + r - q}); ++t)
						d += rf_half(k, r - q) *
							 (powm1(n - l - k) * (1L << (r - 2 * q)))
							 /
							 (Expression(factorial(q)) * factorial(r - 2 * q) *
							  (m + 1 + r))
							 * binomial(Integer(r - q), t)
							 * binomial(Integer(k - 1), m + r)
							 * binomial(Integer(m + 1), 2 * p + 1)
							 * binomial_rational(add(rational(k, 2),
													 integer(r - q + l - k - t -
															 1)), l - k - t)
							 * binomial(Integer(p),
										n - (m + r - q + l - k - t - p));

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

rc d_phi2(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d = Expression(0);
	for (int r = 0; r <= k - 1; ++r)
		for (int m = 0; m <= k - 1 - r; ++m)
			for (int q = 0; q <= r / 2; ++q)
				for (int p = 0; p <= m / 2; ++p)
					for (int t = l - n - k + m + r - p - q;
						 t <= std::min(l - k, r - q); ++t)
						d += rf_half(k, r - q) *
							 (powm1(n - l - k) * (1L << (r - 2 * q)))
							 /
							 (Expression(factorial(q)) * factorial(r - 2 * q) *
							  (m + 1 + r))
							 * binomial(Integer(r - q), t)
							 * binomial(Integer(k - 1), m + r)
							 * binomial(Integer(m + 1), 2 * p + 1)
							 * binomial_rational(add(rational(k, 2),
													 integer(r - q + l - k - t -
															 1)), l - k - t)
							 *
							 binomial_rational(add(rational(1, 2), integer(p)),
											   n - (m + r - q + l - k - t - p));

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

rc c_phi(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d(0);
	for (int r = 0; r <= k - 1; ++r)
		for (int m = 0; m <= k - 1 - r; ++m)
			for (int p = 0; p <= m / 2; ++p)
				for (int q = 0; q <= r / 2; ++q)
					for (int t = 0; t <= std::min(l - k, r - q); ++t) {
						const int w = m + r - q + l - k - t - p;
						Expression h(0);
						for (int i = std::max(0, p - n + 1);
							 i <= std::min(p, p - n + 1 + w); ++i) {
							const int j = w + p - i + 1 - n;
							h += Expression(
									binomial(SymEngine::Integer(2 * p + 1), i))
								 * binomial(SymEngine::Integer(1 + 2 * w), j)
								 * powm1(j);
						}
						for (int i = p - w + n; i <= p; ++i) {
							const int j = w - p + i - n;
							h += Expression(
									binomial(SymEngine::Integer(2 * p + 1), i))
								 * binomial(SymEngine::Integer(1 + 2 * w), j)
								 * powm1(j);
						}
						for (int i = p - w - n; i <= p - n; ++i) {
							const int j = w - p + i + n;
							h -= Expression(
									binomial(SymEngine::Integer(2 * p + 1), i))
								 * binomial(SymEngine::Integer(1 + 2 * w), j)
								 * powm1(j);
						}
						d += h * rf_half(k, r - q) * powm1(l - k)
							 / (Expression(SymEngine::factorial(q)) *
								SymEngine::factorial(r - 2 * q) *
								Expression(m + 1 + r) *
								(1u << (2 * (m + l - k - t) + r + 1)))
							 * binomial(Integer(r - q), t)
							 * binomial(Integer(k - 1), m + r)
							 * binomial(Integer(m + 1), 2 * p + 1)
							 * binomial_rational(add(rational(k, 2),
													 integer(r - q + l - k - t -
															 1)), l - k - t);
					}

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

Expression d_phi_pow_polynomial(int n, int k, int i) {
	assert(n >= 0 && k >= 0 && i >= 0);
	static auto cache = UintsCache<Expression>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) i))
		return *ret;

	// Series for b_{n,i}.
	std::shared_ptr<SeriesBase> series = double_series_power_coeff(n, i);

	// Take the k-th coefficient.
	Expression b_ni_k = series->getItem(k);

	// Substitute a_{n,k} with its inner series using d_phi and n_offset = 1.
	auto d = a_nk_sub(b_ni_k, 1, d_phi);

	cache.insert(d, (uint) n, (uint) k, (uint) i);
	return d;
}

rc d_phi_pow(int n, int k, int l, int i) {
	assert(n >= 0 && k >= 0 && l >= 0 && i >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l, (uint) i))
		return *ret;

	Expression poly = d_phi_pow_polynomial(n, k, i);
	auto d = coeff_of(expand(poly).get_basic(), e2sym, l);

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l, (uint) i);
	return ret;
}

rc d_sin(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d(0);
	const int i_max = std::min(k, 2 * n + 1);
	for (int i = 1; i <= i_max; ++i) {
		const int floor_i_2 = i / 2;
		const int ceil_i_2 = (i + 1) / 2;
		for (int j = std::max(0, ceil_i_2 - l + n);
			 j <= std::min(ceil_i_2, n - floor_i_2); ++j)
			d += Expression(powm1(floor_i_2 + j))
				 * binomial(Integer(ceil_i_2), j)
				 * rc_expr(d_phi_pow(n - floor_i_2 - j, k, l, i))
				 / SymEngine::factorial(i);
	}

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

/** Polynomial for the δ^k coefficient in b_{n,i}, using d_sin for inner coefficients. */
Expression d_sin_pow_polynomial(int n, int k, int i) {
	assert(n >= 0 && k >= 0 && i >= 0);
	static auto cache = UintsCache<Expression>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) i))
		return *ret;

	// Series for b_{n,i} in terms of {a_0,...,a_n}, then take the k-th coefficient.
	std::shared_ptr<SeriesBase> series = double_series_power_coeff(n, i);
	Expression b_ni_k = series->getItem(k);

	// Substitute a_{n,k} with \sum_{l=max(k,n+0)}^{n+k} d_sin(n,k,l) * e2^l  (offset = 0).
	Expression d(a_nk_sub(b_ni_k, 0, d_sin));

	cache.insert(d, (uint) n, (uint) k, (uint) i);
	return d;
}

rc d_sin_pow(int n, int k, int l, int i) {
	assert(n >= 0 && k >= 0 && l >= 0 && i >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l, (uint) i))
		return *ret;

	Expression poly = d_sin_pow_polynomial(n, k, i);
	Expression d(coeff_of(expand(poly), e2sym, l));

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l, (uint) i);
	return ret;
}

rc c_sin(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d(sin_pow_to_cos_mul(n, k, l, 0, 0, d_sin));

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

rc d_cos(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d(0);
	const int i_max = std::min(k, 2 * n + 1);
	for (int i = 1; i <= i_max; ++i) {
		const int floor_i_2 = i / 2;
		const int ceil_i_2 = (i + 1) / 2;
		for (int j = std::max(0, floor_i_2 - l + n);
			 j <= std::min(floor_i_2, n - ceil_i_2); ++j)
			d += Expression(powm1(ceil_i_2 + j))
				 * binomial(Integer(floor_i_2), j)
				 * rc_expr(d_phi_pow(n - ceil_i_2 - j, k, l, i))
				 / SymEngine::factorial(i);
	}

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

rc c_cos(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d(sin_pow_to_cos_mul(n, k, l, 0, -1, d_cos));

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

rc d_N_nkl(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d(0);
	for (int i = 1; i <= std::min(n, l); ++i)
		for (int j = 0; j <= std::min(2 * i, k); ++j)
			d += binomial_rational(rational(1, 2), i)
				 * binomial(Integer(2 * i), j)
				 * Expression(powm1(i))
				 * rc_expr(d_sin_pow(n - i, k, l - i, j));

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

rc bp_nkl(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d(0);
	for (int i = 1; i <= std::min(n, k / 2); ++i)
		for (int j = 0; j <= std::min(i, n - i); ++j)
			d += Expression(powm1(i + j))
				 * binomial(Integer(i), j)
				 / factorial(2 * i)
				 * rc_expr(d_phi_pow(n - i - j, k, l, 2 * i));

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

rc d_h(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d;
	if (k == 0) {
		d = -rc_expr(d_N_nkl(n, k, l));
	} else {
		d = rc_expr(bp_nkl(n, k + 1, l)) - rc_expr(d_N_nkl(n, k, l));
	}
	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

rc c_h(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);
	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	Expression d(sin_pow_to_cos_mul(n, k, l, 1, 0, d_h));

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

}  // namespace point_to_ellipse_series