#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <gmpxx.h>
#include "lexpr.hpp"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// var_fn that maps a_{j,m} -> mpq_class(j * 10 + m) for easy reasoning.
static mpq_class idx_fn(int j, int m) { return mpq_class(j * 10 + m); }

// ---------------------------------------------------------------------------
// Construction and nullptr semantics
// ---------------------------------------------------------------------------

TEST_CASE("lexpr_var returns non-null", "[lexpr]") {
	LExpr x = lexpr_var(1, 0);
	REQUIRE(x != nullptr);
}

TEST_CASE("lexpr_const returns non-null", "[lexpr]") {
	LExpr c = lexpr_const(mpq_class(3, 2));
	REQUIRE(c != nullptr);
}

TEST_CASE("default LExpr is nullptr (zero)", "[lexpr]") {
	LExpr z{};
	REQUIRE(z == nullptr);
}

// ---------------------------------------------------------------------------
// Arithmetic short-circuit rules
// ---------------------------------------------------------------------------

TEST_CASE("nullptr + x == x", "[lexpr]") {
	LExpr z{};
	LExpr x = lexpr_var(1, 0);
	CHECK((z + x) == x);
	CHECK((x + z) == x);
}

TEST_CASE("nullptr + nullptr == nullptr", "[lexpr]") {
	LExpr z{};
	CHECK((z + z) == nullptr);
}

TEST_CASE("nullptr * x == nullptr", "[lexpr]") {
	LExpr z{};
	LExpr x = lexpr_var(1, 0);
	CHECK((z * x) == nullptr);
	CHECK((x * z) == nullptr);
}

TEST_CASE("x * zero scalar == nullptr", "[lexpr]") {
	LExpr x = lexpr_var(1, 0);
	CHECK((x * mpq_class(0)) == nullptr);
}

TEST_CASE("non-null + non-null is non-null", "[lexpr]") {
	LExpr a = lexpr_var(1, 0);
	LExpr b = lexpr_var(0, 1);
	REQUIRE((a + b) != nullptr);
}

TEST_CASE("non-null * non-null is non-null", "[lexpr]") {
	LExpr a = lexpr_var(1, 0);
	LExpr b = lexpr_var(0, 1);
	REQUIRE((a * b) != nullptr);
}

// ---------------------------------------------------------------------------
// Evaluation — basic correctness
// ---------------------------------------------------------------------------

TEST_CASE("eval of nullptr returns 0", "[lexpr]") {
	LExpr z{};
	mpq_class r = lexpr_eval(z, idx_fn);
	CHECK(r == 0);
}

TEST_CASE("eval of const returns that constant", "[lexpr]") {
	mpq_class c(7, 3);
	LExpr e = lexpr_const(c);
	CHECK(lexpr_eval(e, idx_fn) == c);
}

TEST_CASE("eval of var calls var_fn", "[lexpr]") {
	LExpr x = lexpr_var(2, 5);
	// idx_fn(2,5) = 2*10+5 = 25
	CHECK(lexpr_eval(x, idx_fn) == mpq_class(25));
}

TEST_CASE("eval of sum: 2*a_{1,0} + 3*a_{0,1}", "[lexpr]") {
	// With idx_fn: a_{1,0}=10, a_{0,1}=1 => 2*10 + 3*1 = 23
	LExpr a = lexpr_var(1, 0) * mpq_class(2);
	LExpr b = lexpr_var(0, 1) * mpq_class(3);
	LExpr expr = a + b;
	CHECK(lexpr_eval(expr, idx_fn) == mpq_class(23));
}

TEST_CASE("eval of product: a_{1,0} * a_{0,1}", "[lexpr]") {
	// idx_fn: a_{1,0}=10, a_{0,1}=1 => 10*1 = 10
	LExpr expr = lexpr_var(1, 0) * lexpr_var(0, 1);
	CHECK(lexpr_eval(expr, idx_fn) == mpq_class(10));
}

TEST_CASE("eval of scalar * var: (3/2) * a_{1,1}", "[lexpr]") {
	// idx_fn: a_{1,1}=11 => 3/2 * 11 = 33/2
	LExpr expr = mpq_class(3, 2) * lexpr_var(1, 1);
	CHECK(lexpr_eval(expr, idx_fn) == mpq_class(33, 2));
}

TEST_CASE("eval respects rational arithmetic", "[lexpr]") {
	// (1/3)*a_{0,0} + (1/6)*a_{0,0} = (1/2)*a_{0,0}
	// idx_fn: a_{0,0}=0 => result = 0
	// Use a_{1,0} = 10 instead: (1/3)*10 + (1/6)*10 = 10/3 + 10/6 = 30/6 = 5
	LExpr v = lexpr_var(1, 0);
	LExpr expr = (mpq_class(1, 3) * v) + (mpq_class(1, 6) * v);
	CHECK(lexpr_eval(expr, idx_fn) == mpq_class(5));
}

// ---------------------------------------------------------------------------
// DAG memoization: shared node evaluated exactly once
// ---------------------------------------------------------------------------

TEST_CASE("memoization: shared node evaluated once", "[lexpr]") {
	int call_count = 0;
	auto counting_fn = [&](int j, int m) -> mpq_class {
		++call_count;
		return mpq_class(j + m);
	};

	LExpr x = lexpr_var(1, 0);
	// x + x shares the same LNode pointer.
	LExpr expr = x + x;

	std::unordered_map<const LNode*, mpq_class> memo;
	lexpr_eval(expr, counting_fn, memo);

	// var_fn should be called exactly once despite x appearing twice.
	CHECK(call_count == 1);
}

TEST_CASE("memoization: independent vars each evaluated once", "[lexpr]") {
	int call_count = 0;
	auto counting_fn = [&](int j, int m) -> mpq_class {
		++call_count;
		return mpq_class(j + m);
	};

	LExpr x = lexpr_var(1, 0);
	LExpr y = lexpr_var(0, 1);
	LExpr expr = x + y;

	std::unordered_map<const LNode*, mpq_class> memo;
	lexpr_eval(expr, counting_fn, memo);

	CHECK(call_count == 2);
}

TEST_CASE("operator+= works", "[lexpr]") {
	LExpr acc{};
	acc += lexpr_var(1, 0);
	acc += lexpr_var(0, 1);
	// idx_fn: 10 + 1 = 11
	CHECK(lexpr_eval(acc, idx_fn) == mpq_class(11));
}
