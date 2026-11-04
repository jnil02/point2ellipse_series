#pragma once

/*
 * Caching utilities for function taking unsigned integer arguments.
 */

#include <concepts>       // std::unsigned_integral
#include <cstdint>        // std::uint64_t
#include <stdexcept>      // std::overflow_error
#include <unordered_map>  // std::unordered_map

// Cantor pairing, i.e. bijection between N2 and N.
//
// Throws std::overflow_error if the key does not fit in 64 bits; a wrapped key
// would make the cache silently return the value of another argument tuple.
// Chained pairing grows roughly doubly exponentially with the number of
// arguments: with all arguments <= N, the key fits for N <= 234 with four
// arguments and N <= 55108 with three. The check is always on (not only in
// debug builds), as overflow is reached in long high-order (release) runs and
// costs little next to the hash lookup and rational arithmetic.
constexpr std::uint64_t cantor_pairing_two(std::uint64_t a, std::uint64_t b) {
	std::uint64_t s, t, key;
	if (__builtin_add_overflow(a, b, &s)
		|| s == UINT64_MAX  // s + 1 would wrap to 0.
		// t = s(s+1)/2, halving the even factor first so that only the
		// triangular number itself (not s(s+1)) must fit in 64 bits.
		|| (s % 2 == 0 ? __builtin_mul_overflow(s / 2, s + 1, &t)
		               : __builtin_mul_overflow(s, (s + 1) / 2, &t))
		|| __builtin_add_overflow(t, b, &key))
		throw std::overflow_error("cantor_pairing: cache key overflows 64 bits");
	return key;
}

// Chained Cantor pairing for an arbitrary number of unsigned integers.
template<std::unsigned_integral... UInts>
constexpr std::uint64_t cantor_pairing(UInts... args) {
	std::uint64_t arr[] = {static_cast<std::uint64_t>(args)...};
	std::uint64_t key = arr[0];
	for (std::size_t i = 1; i < sizeof...(args); ++i)
		key = cantor_pairing_two(key, arr[i]);
	return key;
}

// Cache for unsigned integers arguments.
template<typename R>
class UintsCache {
public:
	// Get from cache or nullptr.
	template<std::unsigned_integral... UInts>
	R *get(UInts... args) {
		auto key = cantor_pairing(args...);
		auto it = cache.find(key);
		if (it != cache.end())
			return &it->second;
		return nullptr;
	}

	// Insert into cache and return a reference to the stored value.
	template<std::unsigned_integral... UInts>
	const R& insert(R val, UInts... args) {
		auto key = cantor_pairing(args...);
		auto [it, _] = cache.emplace(key, std::move(val));
		return it->second;
	}

private:
	std::unordered_map<std::uint64_t, R> cache;
};
