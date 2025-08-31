
#include <iostream>
#include <iomanip>
#include <mpreal.h>
#include <chrono>

#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/real_mpfr.h>

#include "fourier_series.hpp"
#include "ellipse.hpp"


int main() {
	using mpfr::mpreal;
	using SymEngine::Expression;
	using SymEngine::symbol;
	using SymEngine::real_mpfr;
	using SymEngine::mpfr_class;

	// Set precision for mpreal __before__ using geo constants.
	const mpfr_prec_t bits = 167;  // This is ~50 decimal places.
	set_precision_bits(bits);

//	std::cout.setf(std::ios::scientific);
	std::cout << std::setprecision(50);
	std::cout << std::fixed;              // Force decimal (no scientific).

	// Geodetic coordinates to convert~.
	const mpreal mp_phi = mpreal("43.1") / mpreal("180.0") * mpfr::const_pi();  // Latitude in radians.
	const mpreal mp_h   = mpreal("10000.0");                                    // Altitude in meters.
	const mpreal mp_sin_phi = mpfr::sin(mp_phi);
	const mpreal mp_cos_phi = mpfr::cos(mp_phi);

	// Exact ellipse to Cartesian coordinate transformation.
	const auto [mp_x, mp_y] = mp_ellipse_to_cartesian(mp_phi, mp_h);

	// Derived numeric values.
	const mpreal mp_psi     = mpfr::atan2(mp_y, mp_x);
	const mpreal mp_rho     = mpfr::sqrt(mp_x*mp_x + mp_y*mp_y);
	const mpreal mp_varrho  = mp_a() / mp_rho;
	const mpreal mp_sin_psi = mpfr::sin(mp_psi);
	const mpreal mp_cos_psi = mpfr::cos(mp_psi);

	std::cout << "mp_a: "      << mp_a()    << std::endl;
	std::cout << "mp_phi: "    << mp_phi    << std::endl;
	std::cout << "mp_h: "      << mp_h      << std::endl;
	std::cout << "mp_x: "      << mp_x      << std::endl;
	std::cout << "mp_y: "      << mp_y      << std::endl;
	std::cout << "mp_psi: "    << mp_psi    << std::endl;
	std::cout << "mp_rho: "    << mp_rho    << std::endl;
	std::cout << "mp_varrho: " << mp_varrho << std::endl << std::endl;

	// Cartesian to ellipse reference transformation.
	const auto [mp_phi2, mp_h2] = mp_cartesian_to_ellipse(mp_x, mp_y);

	std::cout << "phi    = " << mp_phi  << std::endl;
	std::cout << "phi2   = " << mp_phi2 << std::endl;
	std::cout << "h      = " << mp_h    << std::endl;
	std::cout << "h2     = " << mp_h2   << std::endl << std::endl;

	// Build and evaluate the series expansions.
	auto start = std::chrono::steady_clock::now();
	const int max_order = 7;

	// Map for substitutions in series.
	// Conversion from mpreal objects via strings.
	SymEngine::map_basic_basic subs;
	subs[symbol("e2")] = real_mpfr(mpfr_class(mp_e2().toString(), bits));
	subs[symbol("varrho")] = real_mpfr(mpfr_class(mp_varrho.toString(), bits));
	subs[symbol("sin_psi")] = real_mpfr(mpfr_class(mp_sin_psi.toString(), bits));
	subs[symbol("cos_psi")] = real_mpfr(mpfr_class(mp_cos_psi.toString(), bits));
	subs[symbol("psi")] = real_mpfr(mpfr_class(mp_psi.toString(), bits));

	std::cout << "phi - psi (ref/pow/mul):\n";
	std::cout << (mp_phi - mp_psi) << "\n";
	std::cout << Expression((phi_in_sin_pow(max_order, max_order) * sin_psi * cos_psi).get_basic()->subs(subs)) << std::endl;
	std::cout << Expression(phi_in_sin_mul(max_order, max_order, max_order).get_basic()->subs(subs)) << std::endl;
	std::cout << Expression((phi_in_sin_pow2(max_order, max_order)).get_basic()->subs(subs)) << std::endl;

	std::cout << "sin(phi)/sin(psi) (ref/pow/mul):\n";
	std::cout << (mp_sin_phi / mp_sin_psi) << "\n";
	std::cout << Expression((sin_phi_in_sin_pow(max_order, max_order) + 1).get_basic()->subs(subs)) << std::endl;
	std::cout << Expression((sin_phi_in_cos_mul(max_order, max_order, max_order) + 1).get_basic()->subs(subs)) << std::endl;
	std::cout << Expression((sin_phi_in_sin_pow2(max_order, max_order, max_order)).get_basic()->subs(subs)) << std::endl;

	std::cout << "cos(phi)/cos(psi) (ref/pow/mul):\n";
	std::cout << (mp_cos_phi / mp_cos_psi) << "\n";
	std::cout << Expression((cos_phi_in_sin_pow(max_order, max_order) + 1).get_basic()->subs(subs)) << std::endl;
	std::cout << Expression((cos_phi_in_cos_mul(max_order, max_order, max_order) + 1).get_basic()->subs(subs)) << std::endl;
	std::cout << Expression((cos_phi_in_sin_pow2(max_order, max_order, max_order)).get_basic()->subs(subs)) << std::endl;

	std::cout << "(h+a-rho)/a (ref/pow/mul):\n";
	std::cout << ((mp_h + mp_a() - mpfr::sqrt(mp_x * mp_x + mp_y * mp_y)) / mp_a()) << std::endl;
	std::cout << Expression((h_in_sin_pow(max_order, max_order)).get_basic()->subs(subs)) << std::endl;
	std::cout << Expression((h_in_cos_mul(max_order, max_order, max_order)).get_basic()->subs(subs)) << std::endl;

	std::cout << "(h+a-rho)/a:\n";
	std::cout << mp_h << std::endl;
	std::cout << mp_rho - mp_a() + mp_a() * rcp_static_cast<const SymEngine::RealMPFR>(Expression((h_in_sin_pow(max_order, max_order)).get_basic()->subs(subs)).get_basic())->as_mpfr().get_mpfr_t() << std::endl;

	auto end = std::chrono::steady_clock::now();     // stop timer
	auto duration = duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Time: " << duration.count() << " ms\n";

	return 0;
}
