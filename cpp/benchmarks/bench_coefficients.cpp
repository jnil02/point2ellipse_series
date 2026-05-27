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
 * Supported functions and their arguments:
 *   ch_evo  n k l       e.g. bench_coefficients ch_evo 10 20 11
 *   cp_evo  n k l       e.g. bench_coefficients cp_evo  8 16  9
 *   a_mr    m r         e.g. bench_coefficients a_mr   12  6
 *   d_phi_evo n k l     e.g. bench_coefficients d_phi_evo 6 3 3
 *   c_phi_evo n k l     e.g. bench_coefficients c_phi_evo 4 8 6
 *
 * Output (one line, tab-separated, parseable):
 *   <function>  <args...>  <wall_ms>
 *
 * Example shell loop for a_mr at increasing orders:
 *   for m in 5 8 10 12 14 16; do
 *       ./bench_coefficients a_mr $m $((m/2))
 *   done
 */

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
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

	if (fn == "d_h_evo") {
		require_args(3);
		auto t0 = Clock::now();
		volatile auto r = d_h_evo(iarg(0), iarg(1), iarg(2));
		(void)r;
		print(elapsed_ms(t0));

	} else if (fn == "cp_evo") {
		require_args(3);
		auto t0 = Clock::now();
		volatile auto r = cp_evo_nkl(iarg(0), iarg(1), iarg(2));
		(void)r;
		print(elapsed_ms(t0));

	} else if (fn == "a_mr") {
		require_args(2);
		auto t0 = Clock::now();
		volatile auto r = a_mr(iarg(0), iarg(1));
		(void)r;
		print(elapsed_ms(t0));

	} else if (fn == "d_phi_evo") {
		require_args(3);
		auto t0 = Clock::now();
		volatile auto r = d_phi_evo(iarg(0), iarg(1), iarg(2));
		(void)r;
		print(elapsed_ms(t0));

	} else if (fn == "c_phi_evo") {
		require_args(3);
		auto t0 = Clock::now();
		volatile auto r = c_phi_evo(iarg(0), iarg(1), iarg(2));
		(void)r;
		print(elapsed_ms(t0));

	} else {
		std::cerr << "Unknown function: " << fn << "\n";
		std::cerr << "Known: c_h_evo, cp_evo, a_mr, d_phi_evo, c_phi_evo\n";
		return 1;
	}

	return 0;
}
