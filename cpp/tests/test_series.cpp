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

#include "expansions_se.hpp"
#include "ellipse/ellipse.hpp"

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
// Fixture: reference values + substitution map for one geodetic test point
// ---------------------------------------------------------------------------

struct Ref {
	mpreal phi, h, x, y, psi, rho, varrho, mp_sin_psi, mp_cos_psi;
	SymEngine::map_basic_basic subs;

	explicit Ref(const std::string& phi_deg) {
		set_precision_bits(BITS);
		std::cout << std::setprecision(50) << std::fixed;

		phi = mpreal(phi_deg) / mpreal("180.0") * mpfr::const_pi();
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
// Far-field series, one geodetic test point per quadrant. phi (the normal
// direction / geodetic latitude) is defined in (-180, 180] deg; |phi| > 90
// places the point at x < 0 (the 2nd / 3rd quadrants), which exercises the
// full-circle Vermeille reference in mp_cartesian_to_ellipse.
// ---------------------------------------------------------------------------

TEST_CASE("far-field series (all quadrants)", "[series]") {
	const std::string phi_deg = GENERATE(std::string("43.1"), std::string("136.9"),
										 std::string("-136.9"), std::string("-43.1"));
	const Ref r(phi_deg);
	INFO("phi = " << phi_deg << " deg");

	// Cartesian <-> ellipse round-trip.
	{
		const auto [phi2, h2] = mp_cartesian_to_ellipse(r.x, r.y);
		assert_close("roundtrip phi [rad]", r.phi, phi2, mpreal("1e-40", BITS));
		assert_close("roundtrip h [m]",     r.h,   h2,   mpreal("1e-40", BITS));
	}

	// phi - psi
	assert_close("phi-psi  sin_pow", r.phi - r.psi,
				 ev(phi_in_sin_pow(MAX_ORD, MAX_ORD) * sin_psi * cos_psi, r.subs), TOL);
	assert_close("phi-psi  sin_mul", r.phi - r.psi,
				 ev(phi_in_sin_mul(MAX_ORD, MAX_ORD, MAX_ORD), r.subs), TOL);
	// phi_in_sin_pow2 absorbs cos(psi) as its positive root sqrt(1-sin^2), so it
	// only represents the right half-plane (cos(psi) >= 0).
	if (r.mp_cos_psi >= 0)
		assert_close("phi-psi  sin_pow2", r.phi - r.psi,
					 ev(phi_in_sin_pow2(MAX_ORD, MAX_ORD), r.subs), mpreal("1e-5", BITS));

	// sin(phi) / sin(psi)
	assert_close("sin(phi)/sin(psi)  sin_pow", mpfr::sin(r.phi) / r.mp_sin_psi,
				 ev(sin_phi_in_sin_pow(MAX_ORD, MAX_ORD), r.subs) + mpreal(1), TOL);
	assert_close("sin(phi)/sin(psi)  cos_mul", mpfr::sin(r.phi) / r.mp_sin_psi,
				 ev(sin_phi_in_cos_mul(MAX_ORD, MAX_ORD, MAX_ORD), r.subs) + mpreal(1), TOL);
	assert_close("sin(phi)/sin(psi)  sin_pow2", mpfr::sin(r.phi) / r.mp_sin_psi,
				 ev(sin_phi_in_sin_pow2(MAX_ORD, MAX_ORD, MAX_ORD), r.subs), TOL);

	// cos(phi) / cos(psi)
	assert_close("cos(phi)/cos(psi)  sin_pow", mpfr::cos(r.phi) / r.mp_cos_psi,
				 ev(cos_phi_in_sin_pow(MAX_ORD, MAX_ORD), r.subs) + mpreal(1), TOL);
	assert_close("cos(phi)/cos(psi)  cos_mul", mpfr::cos(r.phi) / r.mp_cos_psi,
				 ev(cos_phi_in_cos_mul(MAX_ORD, MAX_ORD, MAX_ORD), r.subs) + mpreal(1), TOL);
	assert_close("cos(phi)/cos(psi)  sin_pow2", mpfr::cos(r.phi) / r.mp_cos_psi,
				 ev(cos_phi_in_sin_pow2(MAX_ORD, MAX_ORD, MAX_ORD), r.subs), TOL);

	// (h + a - rho) / a
	assert_close("(h+a-rho)/a  sin_pow", (r.h + mp_a() - r.rho) / mp_a(),
				 ev(h_in_sin_pow(MAX_ORD, MAX_ORD), r.subs), TOL);
	assert_close("(h+a-rho)/a  cos_mul", (r.h + mp_a() - r.rho) / mp_a(),
				 ev(h_in_cos_mul(MAX_ORD, MAX_ORD, MAX_ORD), r.subs), TOL);

	// h in metres (recovered from the h/a series).
	assert_close("h [m]  sin_pow", r.h,
				 ev(h_in_sin_pow(MAX_ORD, MAX_ORD), r.subs) * mp_a() + r.rho - mp_a(), TOL);
}