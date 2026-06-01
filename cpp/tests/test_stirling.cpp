#include <catch.hpp>
#include <gmpxx.h>

#include "stirling.hpp"

using point_to_ellipse_series::stirling1_unsigned;
using point_to_ellipse_series::stirling1_signed;
using point_to_ellipse_series::stirling2;

TEST_CASE("Stirling numbers of the second kind - base cases", "[stirling]")
{
	REQUIRE(stirling2(0, 0) == mpz_class(1));
	REQUIRE(stirling2(5, 0) == mpz_class(0));
	REQUIRE(stirling2(5, 5) == mpz_class(1));
	REQUIRE(stirling2(3, 4) == mpz_class(0));
}

TEST_CASE("Stirling numbers of the second kind - known values", "[stirling]")
{
	REQUIRE(stirling2(5, 1) == mpz_class(1));
	REQUIRE(stirling2(5, 2) == mpz_class(15));
	REQUIRE(stirling2(5, 3) == mpz_class(25));
	REQUIRE(stirling2(5, 4) == mpz_class(10));
}

TEST_CASE("Unsigned Stirling numbers of the first kind - base cases", "[stirling]")
{
	REQUIRE(stirling1_unsigned(0, 0) == mpz_class(1));
	REQUIRE(stirling1_unsigned(5, 0) == mpz_class(0));
	REQUIRE(stirling1_unsigned(5, 5) == mpz_class(1));
	REQUIRE(stirling1_unsigned(3, 4) == mpz_class(0));
}

TEST_CASE("Unsigned Stirling numbers of the first kind - known values", "[stirling]")
{
	REQUIRE(stirling1_unsigned(5, 1) == mpz_class(24));
	REQUIRE(stirling1_unsigned(5, 2) == mpz_class(50));
	REQUIRE(stirling1_unsigned(5, 3) == mpz_class(35));
	REQUIRE(stirling1_unsigned(5, 4) == mpz_class(10));
}

TEST_CASE("Signed Stirling numbers of the first kind - known values", "[stirling]")
{
	REQUIRE(stirling1_signed(5, 1) == mpz_class(24));
	REQUIRE(stirling1_signed(5, 2) == mpz_class(-50));
	REQUIRE(stirling1_signed(5, 3) == mpz_class(35));
	REQUIRE(stirling1_signed(5, 4) == mpz_class(-10));
}

TEST_CASE("Signed/unsigned first kind sign relation", "[stirling]")
{
	for (unsigned n = 1; n <= 8; ++n) {
		for (unsigned k = 1; k <= n; ++k) {
			mpz_class u = stirling1_unsigned(n, k);
			mpz_class s = stirling1_signed(n, k);
			mpz_class expected = ((n - k) % 2 == 0) ? u : -u;
			REQUIRE(s == expected);
		}
	}
}
