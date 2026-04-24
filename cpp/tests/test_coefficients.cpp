#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

#include "coefficients.hpp"

using point_to_ellipse_series::d_phi_evo;
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

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("d_phi_evo matches Python reference", "[coefficients][evo]") {

	// Path is set by CMake via a compile definition TEST_DATA_DIR.
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_phi_evo.csv";

	std::vector<TestRow> rows;
	REQUIRE_NOTHROW(rows = load_csv(csv_path));
	REQUIRE(!rows.empty());

	for (const auto &row : rows) {
		INFO("indices: n=" << row.n << " k=" << row.k << " l=" << row.l);

		rc result = d_phi_evo(row.n, row.k, row.l);
		INFO("c++:    " << row.num << " / " << row.den);
		INFO("python: " << result.num << " / " << result.den);

		// Normalise both fractions before comparing:
		// a/b == c/d  iff  a*d == b*c
		CHECK(result.num * row.den == row.num * result.den);
	}
}

TEST_CASE("d_phi_evo returns zero for invalid indices", "[coefficients][evo]") {
	// l > n/2 + k should always give zero.
	CHECK(d_phi_evo(0, 0, 1).num == 0);
	CHECK(d_phi_evo(2, 1, 3).num == 0);  // n/2+k = 1+1 = 2, l=3 invalid.
	CHECK(d_phi_evo(3, 2, 5).num == 0);  // n/2+k = 1+2 = 3, l=5 invalid.
}
