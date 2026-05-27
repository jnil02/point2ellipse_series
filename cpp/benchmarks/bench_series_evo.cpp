/*
 * Benchmark for inside-evolute series evaluation using mpfr::mpreal.
 *
 * Times a full series evaluation at a concrete numeric point for a given
 * truncation order.  This includes both the coefficient computation (cached
 * within a single run) and the floating-point arithmetic — it is the end-to-end
 * wall time that matters for the convergence sweep.
 *
 * Usage:
 *   bench_series_evo <series> <N> <K> <sin_psi> <rho_ae2> <b_a> [prec_bits]
 *
 * Supported series:
 *   phi_evo_dense      N K sin_psi rho_ae2 b_a
 *   phi_evo_dense_m    M   sin_psi rho_ae2 b_a   (M plays the role of N, K unused)
 *   sin_phi_evo_dense  N K sin_psi rho_ae2 b_a
 *   cos_phi_evo_dense  N K sin_psi rho_ae2 b_a
 *   h_a_evo_dense      N K sin_psi rho_ae2 b_a
 *
 * Optional last argument sets the mpreal working precision in bits (default 256).
 *
 * Output (one line, tab-separated, parseable):
 *   <series>  <N>  <K>  <wall_ms>
 *
 * Example sweep for phi_evo_dense at increasing order:
 *   for N in 2 4 6 8 10 12; do
 *       ./bench_series_evo phi_evo_dense $N $N 0.5 0.3 0.97
 *   done
 */

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>

// mpreal.h must come before fourier_series_evo.hpp, which pulls in series_traits.hpp
// that defines the mpreal specialisations.
#include <mpreal.h>

#include "fourier_series_evo.hpp"

using Clock = std::chrono::steady_clock;

static double elapsed_ms(Clock::time_point t0) {
	return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

int main(int argc, char* argv[]) {
	if (argc < 7) {
		std::cerr << "Usage: bench_series_evo <series> <N> <K> <sin_psi> <rho_ae2> <b_a> [prec_bits]\n";
		return 1;
	}

	const std::string series  = argv[1];
	const int    N        = std::atoi(argv[2]);
	const int    K        = std::atoi(argv[3]);
	const double sin_psi_ = std::atof(argv[4]);
	const double rho_ae2_ = std::atof(argv[5]);
	const double b_a_     = std::atof(argv[6]);
	const int    prec     = (argc >= 8) ? std::atoi(argv[7]) : 256;

	mpfr::mpreal::set_default_prec(prec);

	const mpfr::mpreal sp(sin_psi_);
	const mpfr::mpreal ra(rho_ae2_);
	const mpfr::mpreal ba(b_a_);

	auto print = [&](double ms) {
		std::cout << series << "\t" << N << "\t" << K << "\t" << ms << "\n";
	};

	if (series == "phi_evo_dense") {
		auto t0 = Clock::now();
		volatile auto r = phi_evo_sin_pow_dense<mpfr::mpreal>(N, K, sp, ra, ba);
		(void)r; print(elapsed_ms(t0));

	} else if (series == "phi_evo_dense_m") {
		auto t0 = Clock::now();
		volatile auto r = phi_evo_sin_pow_dense_m<mpfr::mpreal>(N, sp, ra, ba);
		(void)r; print(elapsed_ms(t0));

	} else if (series == "sin_phi_evo_dense") {
		auto t0 = Clock::now();
		volatile auto r = sin_phi_evo_dense<mpfr::mpreal>(N, K, sp, ra, ba);
		(void)r; print(elapsed_ms(t0));

	} else if (series == "cos_phi_evo_dense") {
		auto t0 = Clock::now();
		volatile auto r = cos_phi_evo_dense<mpfr::mpreal>(N, K, sp, ra, ba);
		(void)r; print(elapsed_ms(t0));

	} else if (series == "h_a_evo_dense") {
		auto t0 = Clock::now();
		volatile auto r = h_a_evo_dense<mpfr::mpreal>(N, K, sp, ra, ba);
		(void)r; print(elapsed_ms(t0));

	} else {
		std::cerr << "Unknown series: " << series << "\n";
		std::cerr << "Known: phi_evo_dense, phi_evo_dense_m, sin_phi_evo_dense,\n"
				  << "       cos_phi_evo_dense, h_a_evo_dense\n";
		return 1;
	}

	return 0;
}
