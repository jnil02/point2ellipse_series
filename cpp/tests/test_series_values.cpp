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

TEST_CASE("phi evo series py/cpp", "[series][evo]")
{ check("phi_evo_sin_pow_dense_m2_series.csv", phi_evo_sin_pow_dense_m<double>); }
TEST_CASE("sin(phi) evo series py/cpp", "[series][evo]")
{ check("sin_phi_evo_dense_m_series.csv", sin_phi_evo_dense_m<double>); }
TEST_CASE("cos(phi) evo series py/cpp", "[series][evo]")
{ check("cos_phi_evo_dense_m_series.csv", cos_phi_evo_dense_m<double>); }
TEST_CASE("h evo series py/cpp", "[series][evo]")
{ check("h_a_evo_dense_m2_series.csv", h_a_evo_dense_m3<double>); }