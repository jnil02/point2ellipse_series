
#include <symengine/expression.h>
#include <symengine/ntheory.h>
#include <symengine/rational.h>
#include <symengine/symbol.h>

#include "cache.hpp"
#include "util.hpp"
#include "series_substitution.hpp"

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

using uint = unsigned int;

rc d_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);

	// Only valid (potentially non-zero) for l <= n/2 + k.
	if (l > n / 2 + k)
		return {0, 1};

	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	const int s = n % 2;
	const int m = n + 1 + 2 * k;

	Expression c(0);
	for (int j = 0; j <= l; ++j) {
		const int ka = n / 2 - l + j;

		// Inner sum for b.
		Expression b(0);
		const int i_max = std::min(j, ka);  // Empty range if ka < 0.
		for (int i = 0; i <= i_max; ++i)
			b += Expression(powm1(i))
				 * binomial(Integer(j), (unsigned long) i)
				 * binomial(Integer(ka - i + k), (unsigned long) std::max(0, ka - i));

		// Inner sum for sum_.
		Expression sum_(0);
		for (int q = 2 * j; q <= s + 2 * l; ++q)
			sum_ += Expression(1L << (q - 2 * j))
					* binomial(Integer(q - j), (unsigned long) j)
					* binomial(Integer(n - 1 + 2 * k - q),
							   (unsigned long) (s + 2 * l - q));

		c += sum_ * b;
	}

	// Outer factor: (-1)^(n/2 + l + n + 1) * binomial(m/2, n/2 + k - l) / m
	Expression result =
			Expression(powm1(n / 2 + l + n + 1))
			* binomial_rational(rational(m, 2), (unsigned long) (n / 2 + k - l))
			/ Expression(integer(m))
			* c;

	rc ret = expr_rc(result);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

rc c_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0);

	// Parity constraints — coefficient is zero unless both hold.
	if ((k - n - 1) % 2 != 0 || (l - k) % 2 != 0)
		return {0, 1};

	// Valid range: k >= n+1, 1 <= l <= k, correct parity.
	if (k < n + 1 || l < 1 || l > k)
		return {0, 1};

	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	// Outer binomial factor: binomial(k/2, (k-l)/2) / k.
	// (k-l)/2 >= 0 since (l-k)%2==0 and valid range has l <= k.
	const unsigned long kl_2 = (unsigned long) ((k - l) / 2);
	Expression outer = binomial_rational(rational(k, 2), kl_2)
					   / Expression(integer(k));

	Expression c(0);
	for (int j = 0; j <= (l - 1) / 2; ++j) {

		// Inner sum for b.
		Expression b(0);
		const int i_max = std::min(j, (n + 1 - l + 2 * j) / 2);
		for (int i = 0; i <= i_max; ++i) {
			const int ka = (n - l + 1 + 2 * j - 2 * i) / 2;
			b += Expression(powm1(ka))
				 * binomial(Integer(j), (unsigned long) i)
				 * binomial(Integer(j - i + (k - l) / 2), (unsigned long) ka);
		}

		// Inner sum for s.
		Expression s(0);
		for (int q = 2 * j; q <= l - 1; ++q)
			s += Expression(1L << (q - 2 * j))
				 * Expression(binomial(Integer(k - 2 - q),
									   (unsigned long) (l - 1 - q)))
				 * Expression(binomial(Integer(q - j), (unsigned long) j));

		c += outer * Expression(powm1(l - j)) * s * b;
	}

	rc ret = expr_rc(c);
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

Expression d_phi_pow_evo_polynomial(int n, int k, int i) {
	assert(n >= 0 && k >= 0 && i >= 0);
	static auto cache = UintsCache<Expression>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) i))
		return *ret;

	// Series for b_{n,i}, take the k-th coefficient.
	Expression b_ni_k = double_series_power_coeff(n, i)->getItem(k);

	// Substitute a_{n,k} using a_nk_C with c_phi_evo.
	Expression d = a_nk_sub(b_ni_k, [](int n, int k) {
		return a_nk_C(n, k, c_phi_evo);
	});

	cache.insert(d, (uint) n, (uint) k, (uint) i);
	return d;
}

rc d_phi_pow_evo(int n, int k, int l, int i) {
	assert(n >= 0 && k >= 0 && l >= 0 && i >= 0);

	// Parity constraints from underlying c_phi_evo coefficients.
	if ((i + n - k) % 2 != 0 || (l - k) % 2 != 0)
		return {0, 1};

	static auto cache = UintsCache<rc>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l, (uint) i))
		return *ret;

	Expression poly = d_phi_pow_evo_polynomial(n, k, i);
	auto d = coeff_of(expand(poly).get_basic(), e2sym, l);

	rc ret = expr_rc(d);
	cache.insert(ret, (uint) n, (uint) k, (uint) l, (uint) i);
	return ret;
}

}  // namespace point_to_ellipse_series
