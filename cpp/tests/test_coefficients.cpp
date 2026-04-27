#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

#include <symengine/expression.h>

#include "coefficients.hpp"

using point_to_ellipse_series::d_phi;
using point_to_ellipse_series::d_cos;
using point_to_ellipse_series::d_sin;
using point_to_ellipse_series::d_h;
using point_to_ellipse_series::d_phi_evo;
using point_to_ellipse_series::c_phi_evo;
using point_to_ellipse_series::rc;

// ---------------------------------------------------------------------------
// CSV helpers
// ---------------------------------------------------------------------------

struct TestRow {
	int n, k, l;
	long num, den;
};

static std::vector<TestRow> load_csv(const std::string &path) {
	std::ifstream f(path);
	if (!f.is_open())
		throw std::runtime_error("Could not open test data file: " + path);

	std::vector<TestRow> rows;
	std::string line;

	// Skip header.
	std::getline(f, line);

	while (std::getline(f, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string tok;
		TestRow row{};
		std::getline(ss, tok, ','); row.n   = std::stoi(tok);
		std::getline(ss, tok, ','); row.k   = std::stoi(tok);
		std::getline(ss, tok, ','); row.l   = std::stoi(tok);
		std::getline(ss, tok, ','); row.num = std::stol(tok);
		std::getline(ss, tok, ','); row.den = std::stol(tok);
		rows.push_back(row);
	}
	return rows;
}

static void check_against_csv(const std::string &csv_path,
							  const std::function<rc(int, int, int)> &fn) {
	std::vector<TestRow> rows;
	REQUIRE_NOTHROW(rows = load_csv(csv_path));
	REQUIRE(!rows.empty());

	for (const auto &row : rows) {
		INFO("indices: n=" << row.n << " k=" << row.k << " l=" << row.l);

		rc result = fn(row.n, row.k, row.l);
		INFO("python: " << row.num << " / " << row.den);
		INFO("c++:    " << result.num << " / " << result.den);

		// Normalise both fractions before comparing:
		// a/b == c/d  iff  a*d == b*c
		CHECK(result.num * row.den == row.num * result.den);
	}
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("SymEngine binomial(-1, 0) returns 1", "[symengine]") {
	using SymEngine::Expression;
	using SymEngine::Integer;
	using SymEngine::binomial;
	Expression result = Expression(binomial(Integer(-1), (unsigned long) 0));
	// binomial(n, 0) = 1 for any n by convention.
	CHECK(result == Expression(1));
}

TEST_CASE("d_phi matches Python reference", "[coefficients]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_phi.csv";
	check_against_csv(csv_path, d_phi);
}

TEST_CASE("d_cos matches Python reference", "[coefficients]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_cos.csv";
	check_against_csv(csv_path, d_cos);
}

TEST_CASE("d_sin matches Python reference", "[coefficients]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_sin.csv";
	check_against_csv(csv_path, d_sin);
}

TEST_CASE("d_h matches Python reference", "[coefficients]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_h.csv";
	check_against_csv(csv_path, d_h);
}

TEST_CASE("d_phi_evo matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_phi_evo.csv";
	check_against_csv(csv_path, d_phi_evo);
}

TEST_CASE("c_phi_evo matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/c_phi_evo.csv";
	check_against_csv(csv_path, c_phi_evo);
}
