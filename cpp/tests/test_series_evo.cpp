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

#include "expansions_evo_se.hpp"
#include "ellipse/ellipse.hpp"

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
// Fixture: reference values + substitution map for one polar test point
// ---------------------------------------------------------------------------

struct RefEvo {
	mpreal psi, rho, phi, h, sgn, abs_sin_psi, cos_psi;
	SymEngine::map_basic_basic subs;

	explicit RefEvo(const std::string& psi_deg) {
		set_precision_bits(BITS);
		std::cout << std::setprecision(50) << std::fixed;

		psi = mpreal(psi_deg) / mpreal("180.0") * mpfr::const_pi();
		rho = mpreal("5000.0");

		const mpreal x = rho * mpfr::cos(psi);
		const mpreal y = rho * mpfr::sin(psi);

		auto pair = mp_cartesian_to_ellipse(x, y);
		phi = pair.first;
		h   = pair.second;

		// psi is defined in (-180, 180] deg, so sign(psi) == sign(sin psi).
		sgn         = (psi > 0) ? mpreal(1) : mpreal(-1);
		abs_sin_psi = mpfr::abs(mpfr::sin(psi));
		// Signed cos(psi): its sign carries the x-sign, so the reconstructions
		// below are correct in every quadrant without an explicit reflection.
		cos_psi     = mpfr::cos(psi);

		const mpreal rho_ae2_val = rho / (mp_a() * mp_e2());
		const mpreal b_a_val     = mp_b() / mp_a();

		subs[symbol("sin_psi")] = real_mpfr(mpfr_class(abs_sin_psi.toString(), BITS));
		subs[symbol("rho_ae2")] = real_mpfr(mpfr_class(rho_ae2_val.toString(), BITS));
		subs[symbol("b_a")]     = real_mpfr(mpfr_class(b_a_val.toString(),     BITS));
	}
};

// ---------------------------------------------------------------------------
// Inside-evolute series, one polar test point per quadrant. psi (the polar
// angle) is defined in (-180, 180] deg (see RefEvo), which keeps the quadrant
// sign sgn = sign(sin psi) valid and exercises all four quadrants.
// ---------------------------------------------------------------------------

TEST_CASE("inside-evolute series (all quadrants)", "[series_evo]") {
	const std::string psi_deg = GENERATE(std::string("42.0"), std::string("138.0"),
										 std::string("-138.0"), std::string("-42.0"));
	const RefEvo r(psi_deg);
	INFO("psi = " << psi_deg << " deg");

	// (phi - sgn*pi/2) / (sgn * cos(psi)): the signed cos(psi) makes this
	// correct in every quadrant, so the series compares directly.
	assert_close("phi  evo sparse",
				 (r.phi - r.sgn * mpfr::const_pi() / 2) / (r.sgn * r.cos_psi),
				 ev(phi_evo_sparse(MAX_ORD, MAX_ORD), r.subs), TOL);
	assert_close("phi  evo dense",
				 (r.phi - r.sgn * mpfr::const_pi() / 2) / (r.sgn * r.cos_psi),
				 ev(phi_evo_dense(MAX_ORD), r.subs), TOL);

	// (sin(phi) - sgn) / sgn  (sin(phi) is invariant under the reflection).
	assert_close("(sin(phi)-sgn)/sgn  evo dense", (mpfr::sin(r.phi) - r.sgn) / r.sgn,
				 ev(sin_phi_evo_dense(MAX_ORD), r.subs), TOL);

	// cos(phi): series is cos(phi)/cos(psi); multiply by signed cos(psi).
	assert_close("cos(phi)  evo dense", mpfr::cos(r.phi),
				 r.cos_psi * ev(cos_phi_evo_dense(MAX_ORD), r.subs), TOL);

	// h in metres.
	assert_close("h [m]  evo sparse", r.h,
				 ev(h_evo_sparse(MAX_ORD, MAX_ORD), r.subs) * mp_a() + r.rho * r.abs_sin_psi,
				 TOL * mp_a());
	assert_close("h [m]  evo dense", r.h,
				 ev(h_evo_dense(MAX_ORD), r.subs) * mp_a() - mp_b() + r.rho * r.abs_sin_psi,
				 TOL * mp_a());
}