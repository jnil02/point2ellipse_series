#pragma once

/*
 * Caching utilities for function taking unsigned integer arguments.
 */

#include <concepts>       // std::unsigned_integral
#include <unordered_map>  // std::unordered_map

// Cantor pairing, i.e. bijection between N2 and N.
constexpr long cantor_pairing_two(long a, long b) {
	return ((a + b) * (a + b + 1)) / 2 + b;
}

// Chained Cantor pairing for an arbitrary number of unsigned integers.
template<std::unsigned_integral... UInts>
constexpr long cantor_pairing(UInts... args) {
	long arr[] = {static_cast<long>(args)...};
	long key = arr[0];
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

	// Insert into cache.
	template<std::unsigned_integral... UInts>
	void insert(R val, UInts... args) {
		auto key = cantor_pairing(args...);
		cache.emplace(key, std::move(val));
	}

private:
	std::unordered_map<int, R> cache;
};
