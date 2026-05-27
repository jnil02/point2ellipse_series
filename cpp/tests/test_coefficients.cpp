#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <string>
#include <vector>
#include <stdexcept>

#include "util.hpp"
#include "coefficients.hpp"
#include "coefficients_evo.hpp"

using point_to_ellipse_series::mpz_to_str;
using point_to_ellipse_series::d_phi;
using point_to_ellipse_series::d_cos;
using point_to_ellipse_series::d_sin;
using point_to_ellipse_series::d_h;
using point_to_ellipse_series::d_phi_evo;
using point_to_ellipse_series::d_sin_phi_evo;
using point_to_ellipse_series::d_cos_phi_evo;
using point_to_ellipse_series::c_phi_evo;
using point_to_ellipse_series::c_sin_phi_evo;
using point_to_ellipse_series::c_cos_phi_evo;
using point_to_ellipse_series::c_sin_phi_inv_evo;
using point_to_ellipse_series::a_mr;
using point_to_ellipse_series::B_rt;
using point_to_ellipse_series::C_mt;
using point_to_ellipse_series::R;
using point_to_ellipse_series::B_p;
using point_to_ellipse_series::cp_evo_nkl;
using point_to_ellipse_series::c_h_evo;

// ---------------------------------------------------------------------------
// CSV helpers
// ---------------------------------------------------------------------------

struct TestRow3 {
	int n, k, l;
	long num, den;
};

struct TestRow4 {
	int n, k, l, i;
	long num, den;
};

struct TestRow2 {
	int n, k;
	long num, den;
};

static std::vector<TestRow3> load_csv_3(const std::string &path) {
	std::ifstream f(path);
	if (!f.is_open())
		throw std::runtime_error("Could not open test data file: " + path);

	std::vector<TestRow3> rows;
	std::string line;

	// Skip header.
	std::getline(f, line);

	while (std::getline(f, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string tok;
		TestRow3 row{};
		std::getline(ss, tok, ','); row.n   = std::stoi(tok);
		std::getline(ss, tok, ','); row.k   = std::stoi(tok);
		std::getline(ss, tok, ','); row.l   = std::stoi(tok);
		std::getline(ss, tok, ','); row.num = std::stol(tok);
		std::getline(ss, tok, ','); row.den = std::stol(tok);
		rows.push_back(row);
	}
	return rows;
}

static std::vector<TestRow4> load_csv_4(const std::string &path) {
	std::ifstream f(path);
	if (!f.is_open())
		throw std::runtime_error("Could not open test data file: " + path);

	std::vector<TestRow4> rows;
	std::string line;

	// Skip header.
	std::getline(f, line);

	while (std::getline(f, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string tok;
		TestRow4 row{};
		std::getline(ss, tok, ','); row.n   = std::stoi(tok);
		std::getline(ss, tok, ','); row.k   = std::stoi(tok);
		std::getline(ss, tok, ','); row.l   = std::stoi(tok);
		std::getline(ss, tok, ','); row.i   = std::stoi(tok);
		std::getline(ss, tok, ','); row.num = std::stol(tok);
		std::getline(ss, tok, ','); row.den = std::stol(tok);
		rows.push_back(row);
	}
	return rows;
}

static std::vector<TestRow2> load_csv_2(const std::string &path) {
	std::ifstream f(path);
	if (!f.is_open())
		throw std::runtime_error("Could not open test data file: " + path);

	std::vector<TestRow2> rows;
	std::string line;

	// Skip header.
	std::getline(f, line);

	while (std::getline(f, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string tok;
		TestRow2 row{};
		std::getline(ss, tok, ','); row.n   = std::stoi(tok);
		std::getline(ss, tok, ','); row.k   = std::stoi(tok);
		std::getline(ss, tok, ','); row.num = std::stol(tok);
		std::getline(ss, tok, ','); row.den = std::stol(tok);
		rows.push_back(row);
	}
	return rows;
}

static void check_against_csv_3(const std::string &csv_path,
								const std::function<mpq_class(int, int, int)> &fn) {
	std::vector<TestRow3> rows;
	REQUIRE_NOTHROW(rows = load_csv_3(csv_path));
	REQUIRE(!rows.empty());

	for (const auto &row : rows) {
		INFO("indices: n=" << row.n << " k=" << row.k << " l=" << row.l);
		INFO(csv_path);

		mpq_class result = fn(row.n, row.k, row.l);
		INFO("python: " << row.num << " / " << row.den);
		INFO("c++:    " << mpz_to_str(result.get_num()) << " / " << mpz_to_str(result.get_den()));
		// Normalise both fractions before comparing:
		// a/b == c/d  iff  a*d == b*c
		mpz_class lhs = result.get_num() * row.den;
		mpz_class rhs = mpz_class(row.num) * result.get_den();
		CHECK(lhs == rhs);
	}
}

static void check_against_csv_4(const std::string &csv_path,
								const std::function<mpq_class(int, int, int, int)> &fn) {
	std::vector<TestRow4> rows;
	REQUIRE_NOTHROW(rows = load_csv_4(csv_path));
	REQUIRE(!rows.empty());

	for (const auto &row : rows) {
		INFO("indices: n=" << row.n << " k=" << row.k << " l=" << row.l << " i=" << row.i);
		INFO(csv_path);

		mpq_class result = fn(row.n, row.k, row.l, row.i);
		INFO("python: " << row.num << " / " << row.den);
		INFO("c++:    " << mpz_to_str(result.get_num()) << " / " << mpz_to_str(result.get_den()));

		// Normalise both fractions before comparing:
		// a/b == c/d  iff  a*d == b*c
		mpz_class lhs = result.get_num() * row.den;
		mpz_class rhs = mpz_class(row.num) * result.get_den();
		CHECK(lhs == rhs);
	}
}

static void check_against_csv_2(const std::string &csv_path,
								const std::function<mpq_class(int, int)> &fn) {
	std::vector<TestRow2> rows;
	REQUIRE_NOTHROW(rows = load_csv_2(csv_path));
	REQUIRE(!rows.empty());

	for (const auto &row : rows) {
		INFO("indices: n=" << row.n << " k=" << row.k);
		INFO(csv_path);

		mpq_class result = fn(row.n, row.k);
		INFO("python: " << row.num << " / " << row.den);
		INFO("c++:    " << mpz_to_str(result.get_num()) << " / " << mpz_to_str(result.get_den()));

		// Normalise both fractions before comparing:
		// a/b == c/d  iff  a*d == b*c
		mpz_class lhs = result.get_num() * row.den;
		mpz_class rhs = mpz_class(row.num) * result.get_den();
		CHECK(lhs == rhs);
	}
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("d_phi matches Python reference", "[coefficients]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_phi.csv";
	check_against_csv_3(csv_path, d_phi);
}

TEST_CASE("d_cos matches Python reference", "[coefficients]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_cos.csv";
	check_against_csv_3(csv_path, d_cos);
}

TEST_CASE("d_sin matches Python reference", "[coefficients]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_sin.csv";
	check_against_csv_3(csv_path, d_sin);
}

TEST_CASE("d_h matches Python reference", "[coefficients]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_h.csv";
	check_against_csv_3(csv_path, d_h);
}

TEST_CASE("d_phi_evo matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_phi_evo.csv";
	check_against_csv_3(csv_path, d_phi_evo);
}

TEST_CASE("c_phi_evo matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/c_phi_evo.csv";
	check_against_csv_3(csv_path, c_phi_evo);
}

TEST_CASE("c_sin_phi_evo matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/c_sin_phi_evo.csv";
	check_against_csv_3(csv_path, c_sin_phi_evo);
}

TEST_CASE("d_sin_phi_evo matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_sin_phi_evo.csv";
	check_against_csv_3(csv_path, d_sin_phi_evo);
}

TEST_CASE("c_cos_phi_evo matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/c_cos_phi_evo.csv";
	check_against_csv_3(csv_path, c_cos_phi_evo);
}

TEST_CASE("d_cos_phi_evo matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/d_cos_phi_evo.csv";
	check_against_csv_3(csv_path, d_cos_phi_evo);
}

TEST_CASE("c_sin_phi_inv_evo matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/c_sin_phi_inv_evo.csv";
	check_against_csv_3(csv_path, c_sin_phi_inv_evo);
}

TEST_CASE("a_mr matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/a_mr.csv";
	check_against_csv_2(csv_path, a_mr);
}

TEST_CASE("B_rt matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/B_rt.csv";
	check_against_csv_2(csv_path, B_rt);
}

TEST_CASE("C_mt matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/C_mt.csv";
	check_against_csv_2(csv_path, C_mt);
}

TEST_CASE("R matches Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/R.csv";
	check_against_csv_4(csv_path, R);
}

TEST_CASE("B_p Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/B_p.csv";
	check_against_csv_3(csv_path, B_p);
}

TEST_CASE("cp_evo_nkl Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/cp_evo_nkl.csv";
	check_against_csv_3(csv_path, cp_evo_nkl);
}

TEST_CASE("c_h_evo Python reference", "[coefficients][evo]") {
	const std::string csv_path = std::string(TEST_DATA_DIR) + "/c_h_evo.csv";
	check_against_csv_3(csv_path, c_h_evo);
}
