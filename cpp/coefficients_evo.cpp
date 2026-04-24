
#include <symengine/expression.h>
#include <symengine/ntheory.h>
#include <symengine/rational.h>
#include <symengine/symbol.h>

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

}  // namespace point_to_ellipse_series
