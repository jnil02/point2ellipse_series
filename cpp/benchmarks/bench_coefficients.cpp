/*
 * Micro-benchmark for rational coefficient computation.
 *
 * Measures wall-clock time for one cold call to a named coefficient function.
 * Since each function uses a static in-memory cache, "cold" means the first
 * call in a fresh process.  Run the benchmark repeatedly (from a shell loop)
 * at increasing argument values to build a performance baseline.
 *
 * Usage:
 *   bench_coefficients <function> <args...>
 *
 * Supported functions and their arguments (n, k, l are integers):
 *   d_phi_evo         n k l
 *   d_h_evo           n k l
 *   cp_evo            n k l
 *   a_mr              m r
 *   B_rt              r t
 *   C_mt              m t
 *   B_p               n k p
 *   R                 n k l i
 *
 * Output (one line, tab-separated, parseable):
 *   <function>  <args...>  <wall_ms>
 *
 * Example shell loop for d_h_evo at increasing orders:
 *   for n in 4 6 8 10 12; do
 *       ./bench_coefficients d_h_evo $n $((n*2)) $((n+1))
 *   done
 */

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>

#include "coefficients_evo.hpp"

using namespace point_to_ellipse_series;
using Clock = std::chrono::steady_clock;

static double elapsed_ms(Clock::time_point t0) {
	return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: bench_coefficients <function> <args...>\n";
		return 1;
	}

	const std::string fn = argv[1];

	auto require_args = [&](int needed) {
		if (argc < 2 + needed) {
			std::cerr << fn << " requires " << needed << " integer argument(s)\n";
			std::exit(1);
		}
	};

	auto iarg = [&](int i) { return std::atoi(argv[2 + i]); };

	auto print = [&](double ms) {
		std::cout << fn;
		for (int i = 2; i < argc; ++i) std::cout << "\t" << argv[i];
		std::cout << "\t" << ms << "\n";
	};

	// --- d_xxx_evo coefficients ---

	if (fn == "d_phi_evo") {
		require_args(3);
		auto t0 = Clock::now();
		volatile auto r = d_phi_evo(iarg(0), iarg(1), iarg(2));
		(void)r; print(elapsed_ms(t0));

	} else if (fn == "d_h_evo") {
		require_args(3);
		auto t0 = Clock::now();
		volatile auto r = d_h_evo(iarg(0), iarg(1), iarg(2));
		(void)r; print(elapsed_ms(t0));

		// --- intermediate / helper coefficients ---

	} else if (fn == "cp_evo") {
		require_args(3);
		auto t0 = Clock::now();
		volatile auto r = cp_evo_nkl(iarg(0), iarg(1), iarg(2));
		(void)r; print(elapsed_ms(t0));

	} else if (fn == "a_mr") {
		require_args(2);
		auto t0 = Clock::now();
		volatile auto r = a_mr(iarg(0), iarg(1));
		(void)r; print(elapsed_ms(t0));

	} else if (fn == "B_rt") {
		require_args(2);
		auto t0 = Clock::now();
		volatile auto r = B_rt(iarg(0), iarg(1));
		(void)r; print(elapsed_ms(t0));

	} else if (fn == "C_mt") {
		require_args(2);
		auto t0 = Clock::now();
		volatile auto r = C_mt(iarg(0), iarg(1));
		(void)r; print(elapsed_ms(t0));

	} else if (fn == "B_p") {
		require_args(3);
		auto t0 = Clock::now();
		volatile auto r = B_p(iarg(0), iarg(1), iarg(2));
		(void)r; print(elapsed_ms(t0));

	} else if (fn == "R") {
		require_args(4);
		auto t0 = Clock::now();
		volatile auto r = R(iarg(0), iarg(1), iarg(2), iarg(3));
		(void)r; print(elapsed_ms(t0));

	} else {
		std::cerr << "Unknown function: " << fn << "\n";
		std::cerr << "Known: d_phi_evo, d_h_evo, cp_evo, a_mr, B_rt, C_mt, B_p, R\n";
		return 1;
	}

	return 0;
}
