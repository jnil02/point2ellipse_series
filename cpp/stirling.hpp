#pragma once

#include <symengine/integer.h>
#include <symengine/number.h>
#include <symengine/add.h>
#include <symengine/mul.h>

#include "cache.hpp"

namespace point_to_ellipse_series {

using SymEngine::Integer;
using SymEngine::RCP;
using SymEngine::integer;
using SymEngine::rcp_static_cast;

inline UintsCache<RCP<const Integer>> stirling2_cache;
inline UintsCache<RCP<const Integer>> stirling1_unsigned_cache;

inline RCP<const Integer> as_integer(const SymEngine::RCP<const SymEngine::Number> &x)
{
	return rcp_static_cast<const Integer>(x);
}

// Stirling numbers of the second kind: S(n,k)
inline RCP<const Integer> stirling2(unsigned n, unsigned k)
{
	if (k > n) return integer(0);
	if (n == 0 && k == 0) return integer(1);
	if (k == 0) return integer(0);
	if (k == n) return integer(1);

	if (auto *cached = stirling2_cache.get(n, k))
		return *cached;

	auto result_num = SymEngine::addnum(
			SymEngine::mulnum(integer(k), stirling2(n - 1, k)),
			stirling2(n - 1, k - 1)
	);

	auto result = as_integer(result_num);
	stirling2_cache.insert(result, n, k);
	return result;
}

// Unsigned Stirling numbers of the first kind: c(n,k)
inline RCP<const Integer> stirling1_unsigned(unsigned n, unsigned k)
{
	if (k > n) return integer(0);
	if (n == 0 && k == 0) return integer(1);
	if (k == 0) return integer(0);
	if (k == n) return integer(1);

	if (auto *cached = stirling1_unsigned_cache.get(n, k))
		return *cached;

	auto result_num = SymEngine::addnum(
			SymEngine::mulnum(integer(n - 1), stirling1_unsigned(n - 1, k)),
			stirling1_unsigned(n - 1, k - 1)
	);

	auto result = as_integer(result_num);
	stirling1_unsigned_cache.insert(result, n, k);
	return result;
}

// Signed Stirling numbers of the first kind: s(n,k)
inline RCP<const Integer> stirling1_signed(unsigned n, unsigned k)
{
	auto val = stirling1_unsigned(n, k);

	if ((n - k) % 2 == 0)
		return val;

	return as_integer(SymEngine::mulnum(integer(-1), val));
}

} // namespace point_to_ellipse_series