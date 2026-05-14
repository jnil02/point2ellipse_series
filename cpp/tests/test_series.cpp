// #define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

// This is just to get the tests to run in declaration such that the results
// can easily be compared with the Python side.
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

#include "fourier_series.hpp"
#include "ellipse.hpp"

using mpfr::mpreal;
using SymEngine::Expression;
using SymEngine::symbol;
using SymEngine::real_mpfr;
using SymEngine::mpfr_class;

static const mpfr_prec_t BITS    = 167;   // ~50 decimal places
static const int         MAX_ORD = 7;
static const mpreal      TOL("1e-10", BITS);

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

struct Ref {
	mpreal phi, h, x, y, psi, rho, varrho, mp_sin_psi, mp_cos_psi;
	SymEngine::map_basic_basic subs;

	Ref() {
		set_precision_bits(BITS);
		std::cout << std::setprecision(50) << std::fixed;

		phi = mpreal("43.1") / mpreal("180.0") * mpfr::const_pi();
		h   = mpreal("10000.0");

		auto pair = mp_ellipse_to_cartesian(phi, h);
		x = pair.first;
		y = pair.second;

		psi        = mpfr::atan2(y, x);
		rho        = mpfr::sqrt(x*x + y*y);
		varrho     = mp_a() / rho;
		mp_sin_psi = mpfr::sin(psi);
		mp_cos_psi = mpfr::cos(psi);

		subs[symbol("e2")]      = real_mpfr(mpfr_class(mp_e2().toString(),      BITS));
		subs[symbol("varrho")]  = real_mpfr(mpfr_class(varrho.toString(),       BITS));
		subs[symbol("sin_psi")] = real_mpfr(mpfr_class(mp_sin_psi.toString(),   BITS));
		subs[symbol("cos_psi")] = real_mpfr(mpfr_class(mp_cos_psi.toString(),   BITS));
		subs[symbol("psi")]     = real_mpfr(mpfr_class(psi.toString(),          BITS));
	}
};

// ---------------------------------------------------------------------------
// Round-trip
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(Ref, "roundtrip cartesian to ellipse", "[series]") {
	const auto [phi2, h2] = mp_cartesian_to_ellipse(x, y);
	assert_close("roundtrip phi [rad]", phi, phi2, mpreal("1e-40", BITS));
	assert_close("roundtrip h [m]",     h,   h2,   mpreal("1e-40", BITS));
}

// ---------------------------------------------------------------------------
// phi - psi
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(Ref, "phi minus psi sin_pow", "[series]") {
	const mpreal expected = phi - psi;
	const mpreal result   = ev(phi_in_sin_pow(MAX_ORD, MAX_ORD) * sin_psi * cos_psi, subs);
	assert_close("phi-psi  sin_pow", expected, result, TOL);
}

TEST_CASE_METHOD(Ref, "phi minus psi sin_mul", "[series]") {
	const mpreal expected = phi - psi;
	const mpreal result   = ev(phi_in_sin_mul(MAX_ORD, MAX_ORD, MAX_ORD), subs);
	assert_close("phi-psi  sin_mul", expected, result, TOL);
}

TEST_CASE_METHOD(Ref, "phi minus psi sin_pow2", "[series]") {
	const mpreal expected = phi - psi;
	const mpreal result   = ev(phi_in_sin_pow2(MAX_ORD, MAX_ORD), subs);
	// phi_in_sin_pow2 has worse convergence by design; ~2e-7 at order 7.
	assert_close("phi-psi  sin_pow2", expected, result, mpreal("1e-5", BITS));
}

// ---------------------------------------------------------------------------
// sin(phi) / sin(psi)
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(Ref, "sin(phi)/sin(psi) sin_pow", "[series]") {
	const mpreal expected = mpfr::sin(phi) / mp_sin_psi;
	const mpreal result   = ev(sin_phi_in_sin_pow(MAX_ORD, MAX_ORD), subs) + mpreal(1);
	assert_close("sin(phi)/sin(psi)  sin_pow", expected, result, TOL);
}

TEST_CASE_METHOD(Ref, "sin(phi)/sin(psi) cos_mul", "[series]") {
	const mpreal expected = mpfr::sin(phi) / mp_sin_psi;
	const mpreal result   = ev(sin_phi_in_cos_mul(MAX_ORD, MAX_ORD, MAX_ORD), subs) + mpreal(1);
	assert_close("sin(phi)/sin(psi)  cos_mul", expected, result, TOL);
}

TEST_CASE_METHOD(Ref, "sin(phi)/sin(psi) sin_pow2", "[series]") {
	const mpreal expected = mpfr::sin(phi) / mp_sin_psi;
	const mpreal result   = ev(sin_phi_in_sin_pow2(MAX_ORD, MAX_ORD, MAX_ORD), subs);
	assert_close("sin(phi)/sin(psi)  sin_pow2", expected, result, TOL);
}

// ---------------------------------------------------------------------------
// cos(phi) / cos(psi)
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(Ref, "cos(phi)/cos(psi) sin_pow", "[series]") {
	const mpreal expected = mpfr::cos(phi) / mp_cos_psi;
	const mpreal result   = ev(cos_phi_in_sin_pow(MAX_ORD, MAX_ORD), subs) + mpreal(1);
	assert_close("cos(phi)/cos(psi)  sin_pow", expected, result, TOL);
}

TEST_CASE_METHOD(Ref, "cos(phi)/cos(psi) cos_mul", "[series]") {
	const mpreal expected = mpfr::cos(phi) / mp_cos_psi;
	const mpreal result   = ev(cos_phi_in_cos_mul(MAX_ORD, MAX_ORD, MAX_ORD), subs) + mpreal(1);
	assert_close("cos(phi)/cos(psi)  cos_mul", expected, result, TOL);
}

TEST_CASE_METHOD(Ref, "cos(phi)/cos(psi) sin_pow2", "[series]") {
	const mpreal expected = mpfr::cos(phi) / mp_cos_psi;
	const mpreal result   = ev(cos_phi_in_sin_pow2(MAX_ORD, MAX_ORD, MAX_ORD), subs);
	assert_close("cos(phi)/cos(psi)  sin_pow2", expected, result, TOL);
}

// ---------------------------------------------------------------------------
// (h + a - rho) / a
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(Ref, "(h+a-rho)/a sin_pow", "[series]") {
	const mpreal expected = (h + mp_a() - rho) / mp_a();
	const mpreal result   = ev(h_in_sin_pow(MAX_ORD, MAX_ORD), subs);
	assert_close("(h+a-rho)/a  sin_pow", expected, result, TOL);
}

TEST_CASE_METHOD(Ref, "(h+a-rho)/a cos_mul", "[series]") {
	const mpreal expected = (h + mp_a() - rho) / mp_a();
	const mpreal result   = ev(h_in_cos_mul(MAX_ORD, MAX_ORD, MAX_ORD), subs);
	assert_close("(h+a-rho)/a  cos_mul", expected, result, TOL);
}

// ---------------------------------------------------------------------------
// h in metres (recovered from the h/a series)
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(Ref, "h metres sin_pow", "[series]") {
	const mpreal result = ev(h_in_sin_pow(MAX_ORD, MAX_ORD), subs) * mp_a() + rho - mp_a();
	assert_close("h [m]  sin_pow", h, result, TOL);
}