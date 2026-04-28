#include <catch.hpp>

#include "stirling.hpp"

using point_to_ellipse_series::stirling1_unsigned;
using point_to_ellipse_series::stirling1_signed;
using point_to_ellipse_series::stirling2;

TEST_CASE("Stirling numbers of the second kind - base cases", "[stirling]")
{
	REQUIRE(eq(*stirling2(0,0), *SymEngine::integer(1)));
	REQUIRE(eq(*stirling2(5,0), *SymEngine::integer(0)));
	REQUIRE(eq(*stirling2(5,5), *SymEngine::integer(1)));
	REQUIRE(eq(*stirling2(3,4), *SymEngine::integer(0)));
}

TEST_CASE("Stirling numbers of the second kind - known values", "[stirling]")
{
	REQUIRE(eq(*stirling2(5,1), *SymEngine::integer(1)));
	REQUIRE(eq(*stirling2(5,2), *SymEngine::integer(15)));
	REQUIRE(eq(*stirling2(5,3), *SymEngine::integer(25)));
	REQUIRE(eq(*stirling2(5,4), *SymEngine::integer(10)));
}

TEST_CASE("Unsigned Stirling numbers of the first kind - base cases", "[stirling]")
{
	REQUIRE(eq(*stirling1_unsigned(0,0), *SymEngine::integer(1)));
	REQUIRE(eq(*stirling1_unsigned(5,0), *SymEngine::integer(0)));
	REQUIRE(eq(*stirling1_unsigned(5,5), *SymEngine::integer(1)));
	REQUIRE(eq(*stirling1_unsigned(3,4), *SymEngine::integer(0)));
}

TEST_CASE("Unsigned Stirling numbers of the first kind - known values", "[stirling]")
{
	REQUIRE(eq(*stirling1_unsigned(5,1), *SymEngine::integer(24)));
	REQUIRE(eq(*stirling1_unsigned(5,2), *SymEngine::integer(50)));
	REQUIRE(eq(*stirling1_unsigned(5,3), *SymEngine::integer(35)));
	REQUIRE(eq(*stirling1_unsigned(5,4), *SymEngine::integer(10)));
}

TEST_CASE("Signed Stirling numbers of the first kind - known values", "[stirling]")
{
	REQUIRE(eq(*stirling1_signed(5,1), *SymEngine::integer(24)));
	REQUIRE(eq(*stirling1_signed(5,2), *SymEngine::integer(-50)));
	REQUIRE(eq(*stirling1_signed(5,3), *SymEngine::integer(35)));
	REQUIRE(eq(*stirling1_signed(5,4), *SymEngine::integer(-10)));
}

TEST_CASE("Signed/unsigned first kind sign relation", "[stirling]")
{
	for (unsigned n = 1; n <= 8; ++n) {
		for (unsigned k = 1; k <= n; ++k) {

			auto u = stirling1_unsigned(n, k);
			auto s = stirling1_signed(n, k);

			auto expected =
					((n - k) % 2 == 0)
					? SymEngine::rcp_static_cast<const SymEngine::Basic>(u)
					: SymEngine::neg(u);

			REQUIRE(eq(*s, *expected));
		}
	}
}