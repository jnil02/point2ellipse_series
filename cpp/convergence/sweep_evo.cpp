/*
 * Convergence sweep for the inside-evolute series.
 *
 * For each point (psi, rho) on a polar grid the series phi_evo_sin_pow_dense
 * and h_a_evo_dense are evaluated at increasing truncation orders N=K=1..MAX_ORDER
 * using mpreal arithmetic.  The absolute error against the closed-form Vermeille
 * reference is written to a CSV file.  rho intentionally extends beyond the evolute
 * (rho_evo) so that divergence is visible.
 *
 * Ellipse parameters are set via CMake target_compile_definitions
 * (ELLIPSE_A + ELLIPSE_B or ELLIPSE_INV_F).  The default in CMakeLists.txt
 * is the unit ellipse a=1, b=0.5 (b/a=0.5), which has a large evolute convenient
 * for convergence studies.
 *
 * Output columns:
 *   psi_deg   – polar angle in degrees
 *   rho       – radial distance (same units as a)
 *   rho_evo   – evolute radius at this psi (series valid for rho < rho_evo)
 *   N         – truncation order (N = K)
 *   phi_err   – |phi_approx - phi_true|  [radians]
 *   h_err     – |h_approx   - h_true|   [same units as a]
 */

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <mpreal.h>

#include "ellipse.hpp"
#include "fourier_series_accum.hpp"

using mpfr::mpreal;
using mpfr::const_pi;

static constexpr mpfr_prec_t BITS         = 512;// 167;  // ~50 decimal digits; raise together with CSV_SIG_FIGS
static constexpr int         CSV_SIG_FIGS = 15;  // significant figures written to CSV; keep ≥ BITS*log10(2)

/** Evolute radius at polar angle psi.
 *
 * The evolute of x²/a² + y²/b² = 1 satisfies (aX)^(2/3) + (bY)^(2/3) = (a²-b²)^(2/3).
 * For a ray at angle psi, X = rho*|cos(psi)|, Y = rho*|sin(psi)|, so:
 *
 *   rho_evo(psi) = (a²-b²) / [(a|cos(psi)|)^(2/3) + (b|sin(psi)|)^(2/3)]^(3/2)
 *
 * @param psi  Polar angle in radians (0, pi/2).
 * @return     Evolute radius at psi (same units as a).
 */
static mpreal evolute_rho(const mpreal& psi) {
	const mpreal a  = mp_a();
	const mpreal b  = mp_b();
	const mpreal ac = mpfr::abs(a * mpfr::cos(psi));
	const mpreal bs = mpfr::abs(b * mpfr::sin(psi));
	const mpreal sum = mpfr::cbrt(ac * ac) + mpfr::cbrt(bs * bs);
	return (a * a - b * b) / (sum * mpfr::sqrt(sum));
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {
	set_precision_bits(BITS);

	time_t tstart, tend;
	tstart = time(nullptr);

	const int PSI_STEPS = 91;   // PSI_STEPS evenly spaced psi from 0 to 90.
	const int RHO_STEPS = 50;   // RHO_STEPS evenly spaced rho from rho_min to rho_max.
	// MAX_ORDER x PSI_STEP x RHO_STEP old time -> new time (accum).
	// 15x100x100 took 96s -> 16s
	// 16x100x100 took 119s -> 20s
	// 17x100x100 took 136s -> 20s
	// 18x100x100 took 190s -> 23s
	// 19x100x100 took 230s -> 27s
	// 20x100x100 took 272s -> 33s
	// 21x100x100 took ???s -> 44s
	// 22x100x100 took ???s -> 50s
	// 23x100x100 took ???s -> 63s
	// 24x100x100 took ???s -> 75s
	// 25x100x100 took ???s -> 95s
	// 26x100x100 took ???s -> 124s
	// 27x100x100 took ???s -> 173s
	// 28x100x100 took ???s -> 211s
	// 29x100x100 took ???s -> 268s
	// 30x100x100 took ???s -> 377s
	// Only phi
	// 30x90x50 took 22s
	// 35x90x50 took 28s
	// 40x90x50 took 38s
	// 50x90x50 took 71s
	// 70x90x50 took 140s (512 bits tog 342 sekunder)
	// 80x90x50 took 343s
	const int MAX_ORDER = 70;   // N = K = 1 .. MAX_ORDER

	const mpreal a     = mp_a();
	const mpreal b_a_v = mp_b() / a;
	const mpreal ae2   = a * mp_e2();
	const mpreal pi    = const_pi();

	const mpreal rho_max = mpreal("1.5") * a;  // extends beyond the evolute to show divergence
	const mpreal rho_min = mpreal("0.01") * a;

	const std::string fname = std::string(TEST_DATA_DIR) + "/sweep_evo_m.csv";
	std::ofstream out(fname);
	out << "# a=" << a.toString(15) << " b=" << mp_b().toString(15) << "\n";
//	out << "psi_deg,rho,rho_evo,N,phi_err,h_err\n";
	out << "psi_deg,rho,rho_evo,N,phi_err\n";

	for (int i = 0; i < PSI_STEPS; ++i) {
		const mpreal psi_deg = mpreal(90) * mpreal(i) / mpreal(PSI_STEPS-1);
		const mpreal psi     = psi_deg / mpreal(180) * pi;
		const mpreal abs_sin_psi = mpfr::abs(mpfr::sin(psi));
		const mpreal abs_cos_psi = mpfr::abs(mpfr::cos(psi));
		const mpreal sgn         = mpreal(1);  // psi in (0°, 90°) so sin(psi) > 0

		const mpreal rho_evo = evolute_rho(psi);

		std::cerr << "psi = " << psi_deg.toString(4) << "°  rho_evo = " << rho_evo << "\n";

		for (int j = 0; j < RHO_STEPS; ++j) {
			const mpreal rho = rho_min + (rho_max - rho_min) * mpreal(j) / mpreal(RHO_STEPS-1);

			// True solution via Vermeille closed form.
			const mpreal x = rho * mpfr::cos(psi);
			const mpreal y = rho * mpfr::sin(psi);
			auto [true_phi, true_h] = mp_cartesian_to_ellipse(x, y);

			// Shared series inputs.
			const mpreal rho_ae2_v = rho / ae2;

			EvoBasePowers<mpreal> pows(abs_sin_psi, rho_ae2_v, b_a_v, MAX_ORDER);
			PhiEvoAccum<mpreal>   phi_acc(pows);
//			HAEvoAccum<mpreal>    h_acc(pows);

			for (int N = 1; N <= MAX_ORDER; ++N) {
				phi_acc.addOrder(N);
//				h_acc.addOrder(N);

				// phi via phi_evo_sin_pow_dense_m (incremental):
				//   series = (phi - sgn*pi/2) / (sgn*|cos(psi)|)
				//   phi    = sgn*pi/2 + sgn*|cos(psi)| * series
				const mpreal phi_approx = sgn * pi / 2 + sgn * abs_cos_psi * phi_acc.value();
//				const mpreal phi_series = phi_evo_sin_pow_dense_m<mpreal>(
//						N, abs_sin_psi, rho_ae2_v, b_a_v);
//				const mpreal phi_approx = sgn * pi / 2 + sgn * abs_cos_psi * phi_series;


				// h via h_a_evo_dense_m (incremental):
				//   series = h/a  →  h = series * a
//				const mpreal h_approx = h_acc.value() * a;
//				const mpreal h_series = h_a_evo_dense_m<mpreal>(
//						N, abs_sin_psi, rho_ae2_v, b_a_v);
//				const mpreal h_approx = h_series * a;

				// Compute series errors,
				const mpreal phi_err    = mpfr::abs(phi_approx - true_phi);
//				const mpreal h_err    = mpfr::abs(h_approx - true_h);

				out << psi_deg.toString(CSV_SIG_FIGS)  << ","
					<< rho.toString(CSV_SIG_FIGS)      << ","
					<< rho_evo.toString(CSV_SIG_FIGS)  << ","
					<< N                               << ","
					<< phi_err.toString(CSV_SIG_FIGS)
//					<< ","
//					<< h_err.toString(CSV_SIG_FIGS)
					<< "\n";
			}
		}
	}

	std::cout << "Written " << fname << "\n";

	tend = time(nullptr);
	std::cout << "Sweep took "<< difftime(tend, tstart) <<" second(s)."<< std::endl;
	return 0;
}
