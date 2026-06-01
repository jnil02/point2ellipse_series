#pragma once

#include <gmpxx.h>

#include "cache.hpp"

namespace point_to_ellipse_series {

// Stirling numbers of the second kind: S(n,k)
inline mpz_class stirling2(unsigned n, unsigned k)
{
	static UintsCache<mpz_class> cache;
	if (auto *v = cache.get(n, k)) return *v;

	mpz_class result;
	if (n == 0 && k == 0)  result = 1;
	else if (k == 0 || k > n) result = 0;
	else if (k == n)        result = 1;
	else result = mpz_class(k) * stirling2(n - 1, k) + stirling2(n - 1, k - 1);

	cache.insert(result, n, k);
	return result;
}

// Unsigned Stirling numbers of the first kind: c(n,k)
inline mpz_class stirling1_unsigned(unsigned n, unsigned k)
{
	static UintsCache<mpz_class> cache;
	if (auto *v = cache.get(n, k)) return *v;

	mpz_class result;
	if (n == 0 && k == 0)  result = 1;
	else if (k == 0 || k > n) result = 0;
	else if (k == n)        result = 1;
	else result = mpz_class(n - 1) * stirling1_unsigned(n - 1, k) + stirling1_unsigned(n - 1, k - 1);

	cache.insert(result, n, k);
	return result;
}

// Signed Stirling numbers of the first kind: s(n,k)
inline mpz_class stirling1_signed(unsigned n, unsigned k)
{
	mpz_class c = stirling1_unsigned(n, k);
	return ((n - k) % 2 == 0) ? c : -c;
}

} // namespace point_to_ellipse_series
