/*
 * Coefficient ROC diagnostic for the inside-evolute phi series.
 *
 * For a set of fixed ψ angles, computes the per-slab coefficient
 *
 *   C_m(sin ψ, b/a)  =  addOrder(m) evaluated with rho_ae2 = 1
 *
 * i.e. the slab value with the rho_ae2^m factor removed.  By the
 * Cauchy-Hadamard theorem the ratio |C_{m+1}/C_m| converges to 1/R where
 * R is the radius of convergence of the series in rho_ae2.
 *
 *   R = 1/e   → ratio → e   (ae singularity; complex, causes oscillation)
 *   R = ρ_evo/(ae²) → ratio → ae²/ρ_evo  (evolute; real, monotone convergence)
 *
 * For the default ellipse (a=1, b=0.25): e≈0.9682, ae²=0.9375.
 * - ψ < ~62°: evolute is binding (ρ_evo < ae) → ratio < 0.968, monotone
 * - ψ > ~62°: ae singularity is binding         → ratio → 0.968, oscillating
 *
 * Output columns:
 *   psi_deg  – polar angle in degrees
 *   m        – slab order
 *   C_m      – signed slab coefficient
 *   abs_C_m  – |C_m|
 *   ratio    – |C_m| / |C_{m-1}|  (-1 for m=1)
 */

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <mpreal.h>

#include "ellipse.hpp"
#include "fourier_series_accum.hpp"

using mpfr::mpreal;
using mpfr::const_pi;

static constexpr mpfr_prec_t BITS      = 512;
static constexpr int         MAX_ORDER = 200;
static constexpr int         CSV_SIG   = 15;

int main() {
	set_precision_bits(BITS);

	time_t tstart = time(nullptr);

	const mpreal pi   = const_pi();
	const mpreal a    = mp_a();
	const mpreal b_a  = mp_b() / a;

	const std::vector<double> psi_degs_d = {10, 30, 45, 60, 70, 75, 80};

	const std::string fname = std::string(TEST_DATA_DIR) + "/diag_coeff_evo.csv";
	std::ofstream out(fname);
	out << "# a=" << a.toString(15) << " b=" << mp_b().toString(15) << "\n";
	out << "psi_deg,m,C_m,abs_C_m,ratio\n";

	for (double psi_deg_d : psi_degs_d) {
		const mpreal psi     = mpreal(psi_deg_d) / mpreal(180) * pi;
		const mpreal abs_sin = mpfr::abs(mpfr::sin(psi));

		std::cerr << "psi = " << psi_deg_d << "°\n";

		// rho_ae2_v = 1 → pows_.rho_[m] = 1  →  addOrder(m) returns C_m exactly
		EvoBasePowers<mpreal> pows(abs_sin, mpreal(1), b_a, MAX_ORDER);
		PhiEvoAccum<mpreal>   phi_acc(pows);

		mpreal prev_abs(0);
		for (int m = 1; m <= MAX_ORDER; ++m) {
			const mpreal cm     = phi_acc.addOrder(m);
			const mpreal cm_abs = mpfr::abs(cm);
			const mpreal ratio  = (m > 1 && prev_abs > 0)
								  ? cm_abs / prev_abs : mpreal(-1);

			out << psi_deg_d              << ","
				<< m                      << ","
				<< cm.toString(CSV_SIG)   << ","
				<< cm_abs.toString(CSV_SIG) << ","
				<< ratio.toString(CSV_SIG) << "\n";

			prev_abs = cm_abs;
		}
	}

	std::cout << "Written " << fname << "\n";
	time_t tend = time(nullptr);
	std::cout << "Diagnostic took " << difftime(tend, tstart) << " second(s).\n";
	return 0;
}
