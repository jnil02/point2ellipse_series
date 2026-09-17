#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>

#include "fourier_series_evo.hpp"
#include "convergence/fourier_series_accum.hpp"
#include "series_traits/series_traits_double.hpp"

using point_to_ellipse_series::phi_evo_sparse;
using point_to_ellipse_series::phi_evo_dense;
using point_to_ellipse_series::sin_phi_evo_sparse;
using point_to_ellipse_series::sin_phi_evo_dense;
using point_to_ellipse_series::cos_phi_evo_sparse;
using point_to_ellipse_series::cos_phi_evo_dense;
using point_to_ellipse_series::h_evo_sparse;
using point_to_ellipse_series::h_evo_dense;

struct SRow { int order; double sin_psi, rho_ae2, b_a, value; };

static std::vector<SRow> load(const std::string& path) {
	std::ifstream f(path);
	if (!f.is_open())
		throw std::runtime_error("Could not open test data file: " + path);
	std::vector<SRow> rows;
	std::string line;
	std::getline(f, line);                       // skip header
	while (std::getline(f, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string tok;
		SRow r{};
		std::getline(ss, tok, ','); r.order   = std::stoi(tok);
		std::getline(ss, tok, ','); r.sin_psi = std::stod(tok);
		std::getline(ss, tok, ','); r.rho_ae2 = std::stod(tok);
		std::getline(ss, tok, ','); r.b_a     = std::stod(tok);
		std::getline(ss, tok, ','); r.value   = std::stod(tok);
		rows.push_back(r);
	}
	return rows;
}

static void check(const std::string& csv,
				  const std::function<double(int,double,double,double)>& fn) {
	auto rows = load(std::string(TEST_DATA_DIR) + "/" + csv);
	REQUIRE(!rows.empty());
	for (auto& r : rows) {
		double got = fn(r.order, r.sin_psi, r.rho_ae2, r.b_a);
		double tol = 1e-9 * std::max(1.0, std::abs(r.value));   // abs+rel
		INFO(csv << "  py=" << r.value << " cpp=" << got);
		CHECK(std::abs(got - r.value) < tol);
	}
}

static void check(const std::string& csv,
				  const std::function<double(int,int,double,double,double)>& fn) {
	auto rows = load(std::string(TEST_DATA_DIR) + "/" + csv);
	REQUIRE(!rows.empty());
	for (auto& r : rows) {
		double got = fn(r.order, r.order, r.sin_psi, r.rho_ae2, r.b_a);
		double tol = 1e-9 * std::max(1.0, std::abs(r.value));   // abs+rel
		INFO(csv << "  py=" << r.value << " cpp=" << got);
		CHECK(std::abs(got - r.value) < tol);
	}
}

TEST_CASE("phi evo series sparse py/cpp", "[series][evo]")
{ check("phi_evo_sparse_series.csv", phi_evo_sparse<double>); }
TEST_CASE("phi evo series py/cpp", "[series][evo]")
{ check("phi_evo_dense_series.csv", phi_evo_dense<double>); }

TEST_CASE("sin(phi) evo series sparse py/cpp", "[series][evo]")
{ check("sin_phi_evo_sparse_series.csv", sin_phi_evo_sparse<double>); }
TEST_CASE("sin(phi) evo series py/cpp", "[series][evo]")
{ check("sin_phi_evo_dense_series.csv", sin_phi_evo_dense<double>); }

TEST_CASE("cos(phi) evo series sparse py/cpp", "[series][evo]")
{ check("cos_phi_evo_sparse_series.csv", cos_phi_evo_sparse<double>); }
TEST_CASE("cos(phi) evo series py/cpp", "[series][evo]")
{ check("cos_phi_evo_dense_series.csv", cos_phi_evo_dense<double>); }

TEST_CASE("h evo series sparse py/cpp", "[series][evo]")
{ check("h_evo_sparse_series.csv", h_evo_sparse<double>); }
TEST_CASE("h evo series py/cpp", "[series][evo]")
{ check("h_evo_dense_series.csv", h_evo_dense<double>); }

// ---------------------------------------------------------------------------
// Accumulator regression: the incremental accumulators used by sweep_evo /
// diag_coeff_evo must reproduce their dense series exactly.  The CSV `value`
// column is the (Python-verified) dense series value, so summing addOrder over
// all orders must match it.
// ---------------------------------------------------------------------------

TEST_CASE("PhiEvoAccum matches phi_evo_dense", "[series][evo][accum]") {
	auto rows = load(std::string(TEST_DATA_DIR) + "/phi_evo_dense_series.csv");
	REQUIRE(!rows.empty());
	for (auto& r : rows) {
		EvoBasePowers<double> pows(r.sin_psi, r.rho_ae2, r.b_a, r.order);
		PhiEvoAccum<double>   acc(pows);
		for (int N = 1; N <= r.order; ++N) acc.addOrder(N);
		double tol = 1e-9 * std::max(1.0, std::abs(r.value));
		INFO("PhiEvoAccum  series=" << r.value << "  accum=" << acc.value());
		CHECK(std::abs(acc.value() - r.value) < tol);
	}
}

TEST_CASE("HEvoAccum matches h_evo_dense", "[series][evo][accum]") {
	auto rows = load(std::string(TEST_DATA_DIR) + "/h_evo_dense_series.csv");
	REQUIRE(!rows.empty());
	for (auto& r : rows) {
		EvoBasePowers<double> pows(r.sin_psi, r.rho_ae2, r.b_a, r.order);
		HEvoAccum<double>    acc(pows);
		for (int N = 0; N <= r.order; ++N) acc.addOrder(N);
		double tol = 1e-9 * std::max(1.0, std::abs(r.value));
		INFO("HAEvoAccum  series=" << r.value << "  accum=" << acc.value());
		CHECK(std::abs(acc.value() - r.value) < tol);
	}
}