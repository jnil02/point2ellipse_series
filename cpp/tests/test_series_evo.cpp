// #define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

// Run tests in declaration order to match the Python test file.
int main(int argc, char* argv[]) {
	Catch::Session session;
	int result = session.applyCommandLine(argc, argv);
	if (result != 0) return result;
	session.configData().runOrder = Catch::RunTests::InDeclarationOrder;
	return session.run();
}

#include <iostream>
#include <iomanip>
#include <string>

#include <mpreal.h>
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/real_mpfr.h>

#include "fourier_series_evo.hpp"
#include "ellipse.hpp"

using mpfr::mpreal;
using SymEngine::Expression;
using SymEngine::symbol;
using SymEngine::real_mpfr;
using SymEngine::mpfr_class;

static const mpfr_prec_t BITS    = 167;   // ~50 decimal places
static const int         MAX_ORD = 11;
static const mpreal      TOL("1e-9", BITS);

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static mpreal ev(const Expression& expr, const SymEngine::map_basic_basic& s) {
	return mpreal(rcp_static_cast<const SymEngine::RealMPFR>(
			Expression(expr.get_basic()->subs(s)).get_basic()
	)->as_mpfr().get_mpfr_t());
}

static void assert_close(const std::string& title,
						 const mpreal& expected,
						 const mpreal& actual,
						 const mpreal& tol) {
	const mpreal abs_err = abs(actual - expected);
	const mpreal rel_err = (expected != 0) ? abs_err / abs(expected) : mpreal("inf");
	std::cout << "\n" << title << "\n"
			  << "  expected: " << expected << "\n"
			  << "  actual:   " << actual   << "\n"
			  << "  abs err:  " << abs_err  << "  (tol " << tol << ")\n"
			  << "  rel err:  " << rel_err  << "\n";
	REQUIRE(abs_err < tol);
}

// ---------------------------------------------------------------------------
// Fixture: reference values + substitution map
// ---------------------------------------------------------------------------

struct RefEvo {
	mpreal psi, rho, phi, h, sgn, abs_sin_psi, abs_cos_psi;
	SymEngine::map_basic_basic subs;

	RefEvo() {
		set_precision_bits(BITS);
		std::cout << std::setprecision(50) << std::fixed;

		psi = mpreal("138.0") / mpreal("180.0") * mpfr::const_pi();
		rho = mpreal("5000.0");

		const mpreal x = rho * mpfr::cos(psi);
		const mpreal y = rho * mpfr::sin(psi);

		auto pair = mp_cartesian_to_ellipse(x, y);
		phi = pair.first;
		h   = pair.second;

		sgn         = (psi > 0) ? mpreal(1) : mpreal(-1);
		abs_sin_psi = mpfr::abs(mpfr::sin(psi));
		abs_cos_psi = mpfr::abs(mpfr::cos(psi));

		const mpreal rho_ae2_val = rho / (mp_a() * mp_e2());
		const mpreal b_a_val     = mp_b() / mp_a();

		subs[symbol("sin_psi")] = real_mpfr(mpfr_class(abs_sin_psi.toString(), BITS));
		subs[symbol("rho_ae2")] = real_mpfr(mpfr_class(rho_ae2_val.toString(), BITS));
		subs[symbol("b_a")]     = real_mpfr(mpfr_class(b_a_val.toString(),     BITS));
	}
};

// ---------------------------------------------------------------------------
// (phi - sgn * pi/2) / (sgn * |cos(psi)|)
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(RefEvo, "(phi-sgn*pi/2)/(sgn*|cos(psi)|) evo_dense", "[series_evo]") {
	const mpreal expected = (phi - sgn * mpfr::const_pi() / 2) / (sgn * abs_cos_psi);
	const mpreal result   = ev(phi_evo_sin_pow_dense(MAX_ORD, MAX_ORD), subs);
	assert_close("(phi-sgn*pi/2)/(sgn*|cos(psi)|)  evo_dense", expected, result, TOL);
}

TEST_CASE_METHOD(RefEvo, "(phi-sgn*pi/2)/(sgn*|cos(psi)|) evo_dense_m", "[series_evo]") {
	const mpreal expected = (phi - sgn * mpfr::const_pi() / 2) / (sgn * abs_cos_psi);
	const mpreal result   = ev(phi_evo_sin_pow_dense_m(MAX_ORD), subs);
	assert_close("(phi-sgn*pi/2)/(sgn*|cos(psi)|)  evo_dense_m", expected, result, TOL);
}

// ---------------------------------------------------------------------------
// (sin(phi) - sgn) / sgn
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(RefEvo, "(sin(phi)-sgn)/sgn evo_dense", "[series_evo]") {
	const mpreal expected = (mpfr::sin(phi) - sgn) / sgn;
	const mpreal result   = ev(sin_phi_evo_dense(MAX_ORD, MAX_ORD), subs);
	assert_close("(sin(phi)-sgn)/sgn  evo_dense", expected, result, TOL);
}

// ---------------------------------------------------------------------------
// cos(phi) / |cos(psi)|
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(RefEvo, "cos(phi)/|cos(psi)| evo_dense", "[series_evo]") {
	const mpreal expected = mpfr::cos(phi) / abs_cos_psi;
	const mpreal result   = ev(cos_phi_evo_dense(MAX_ORD, MAX_ORD), subs);
	assert_close("cos(phi)/|cos(psi)|  evo_dense", expected, result, TOL);
}

// ---------------------------------------------------------------------------
// h in metres
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(RefEvo, "h metres evo", "[series_evo]") {
	const mpreal result = ev(h_a_evo(MAX_ORD, MAX_ORD), subs) * mp_a();
	assert_close("h [m]  evo", h, result, TOL * mp_a());
}

TEST_CASE_METHOD(RefEvo, "h metres evo_dense", "[series_evo]") {
	const mpreal result = ev(h_a_evo_dense(MAX_ORD, MAX_ORD), subs) * mp_a();
	assert_close("h [m]  evo_dense", h, result, TOL * mp_a());
}
