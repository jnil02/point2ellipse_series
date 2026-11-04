#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "detail/cache.hpp"

#include <cstdint>
#include <set>
#include <stdexcept>

TEST_CASE("cantor_pairing matches the closed form", "[cache]") {
	REQUIRE(cantor_pairing_two(0, 0) == 0);
	REQUIRE(cantor_pairing_two(3, 4) == (7 * 8) / 2 + 4);
	REQUIRE(cantor_pairing(5u) == 5);
}

TEST_CASE("cantor_pairing is injective on a small grid", "[cache]") {
	std::set<std::uint64_t> keys;
	for (unsigned k = 0; k < 20; ++k)
		for (unsigned l = 0; l < 20; ++l)
			for (unsigned n = 0; n < 20; ++n)
				for (unsigned i = 0; i < 20; ++i)
					REQUIRE(keys.insert(cantor_pairing(k, l, n, i)).second);
}

TEST_CASE("cantor_pairing throws instead of wrapping", "[cache]") {
	// Four arguments: the largest equal argument that fits is 234.
	REQUIRE_NOTHROW(cantor_pairing(234u, 234u, 234u, 234u));
	REQUIRE_THROWS_AS(cantor_pairing(235u, 235u, 235u, 235u), std::overflow_error);
	// Three arguments: 55108.
	REQUIRE_NOTHROW(cantor_pairing(55108u, 55108u, 55108u));
	REQUIRE_THROWS_AS(cantor_pairing(55109u, 55109u, 55109u), std::overflow_error);
	// Edge cases in the pairing itself.
	REQUIRE_THROWS_AS(cantor_pairing_two(UINT64_MAX, 0), std::overflow_error);
	REQUIRE_THROWS_AS(cantor_pairing_two(UINT64_MAX / 2, UINT64_MAX / 2 + 1),
					  std::overflow_error);
}

TEST_CASE("UintsCache stores and retrieves by argument tuple", "[cache]") {
	UintsCache<int> cache;
	REQUIRE(cache.get(1u, 2u) == nullptr);
	cache.insert(12, 1u, 2u);
	cache.insert(21, 2u, 1u);
	REQUIRE(*cache.get(1u, 2u) == 12);
	REQUIRE(*cache.get(2u, 1u) == 21);
}
