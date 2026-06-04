#pragma once

/*
 * Incremental (slab-by-slab) accumulators for the inside-evolute series.
 *
 * Evaluating phi_evo_sin_pow_dense_m or h_a_evo_dense_m at each truncation
 * order N in a loop re-sums all prior slabs, giving O(N^2) work per (psi, rho)
 * point and O(N^3) total.  The classes here decouple slab addition so each
 * new order costs only the m=N contribution — O(N^2) total.
 *
 * Usage:
 *
 *   EvoBasePowers<mpreal> pows(sin_psi, rho_ae2, b_a, MAX_ORDER);
 *   PhiEvoAccum<mpreal>   phi_acc(pows);
 *   HAEvoAccum<mpreal>    h_acc(pows);
 *   for (int N = 1; N <= MAX_ORDER; ++N) {
 *       mpreal cm = phi_acc.addOrder(N);  // returns slab C_m; also updates value()
 *       h_acc.addOrder(N);
 *       use(phi_acc.value(), h_acc.value());
 *   }
 *
 * EvoBasePowers holds power tables shared across all accumulators for the same
 * (psi, rho) point; adding new series requires only a new accumulator class.
 */

#include <vector>

#include "fourier_series_evo.hpp"

// ---------------------------------------------------------------------------
// EvoBasePowers — precomputed power tables for the three series bases
// ---------------------------------------------------------------------------

template<typename T>
struct EvoBasePowers {
	std::vector<T> sin_;  // sin_[i] = sin_psi_v^i,  i = 0..max_order
	std::vector<T> rho_;  // rho_[i] = rho_ae2_v^i,  i = 0..max_order
	std::vector<T> b_;    // b_[i]   = b_a_v^i,       i = 0..max_order+1

	EvoBasePowers(const T& sin_psi_v, const T& rho_ae2_v, const T& b_a_v,
				  int max_order) {
		sin_.resize(max_order + 1);
		rho_.resize(max_order + 1);
		b_.resize(max_order + 2);

		sin_[0] = T(1);
		for (int i = 1; i <= max_order; ++i)
			sin_[i] = sin_[i - 1] * sin_psi_v;

		rho_[0] = T(1);
		for (int i = 1; i <= max_order; ++i)
			rho_[i] = rho_[i - 1] * rho_ae2_v;

		b_[0] = T(1);
		for (int i = 1; i <= max_order + 1; ++i)
			b_[i] = b_[i - 1] * b_a_v;
	}
};

// ---------------------------------------------------------------------------
// PhiEvoAccum — incremental accumulator for phi_evo_sin_pow_dense_m
// ---------------------------------------------------------------------------

template<typename T>
class PhiEvoAccum {
public:
	explicit PhiEvoAccum(const EvoBasePowers<T>& pows)
			: pows_(pows), accum_(T(0)) {}

	// Add the m=N slab (mirrors the inner loop of phi_evo_sin_pow_dense_m).
	// Returns the slab value so callers can inspect C_m without a separate class.
	T addOrder(int N) {
		T slab(0);
		const int L      = (N - 1) / 2;
		const int parity = (N - 1) % 2;
		const T&  rho_m  = pows_.rho_[N];
		for (int k = 0; k <= L; ++k) {
			const int n     = N - 1 - 2 * k;
			const T&  sin_n = pows_.sin_[n];
			for (int l = 0; l <= L; ++l)
				slab = slab
					   + point_to_ellipse_series::series_coeff<T>(d_phi_evo(n, k, l))
						 * sin_n * rho_m * pows_.b_[parity + 1 + 2 * l];
		}
		accum_ = accum_ + slab;
		return slab;
	}

	const T& value() const { return accum_; }

private:
	const EvoBasePowers<T>& pows_;
	T accum_;
};

// ---------------------------------------------------------------------------
// HAEvoAccum — incremental accumulator for h_a_evo_dense_m
// ---------------------------------------------------------------------------

template<typename T>
class HAEvoAccum {
public:
	explicit HAEvoAccum(const EvoBasePowers<T>& pows)
			: pows_(pows), accum_(T(0)) {}

	// Add the m=N slab (mirrors the inner loop of h_a_evo_dense_m).
	// Returns the slab value.
	T addOrder(int N) {
		T slab(0);
		const int sn    = N % 2;
		const int lmax  = (N + 1) / 2;
		const T&  rho_m = pows_.rho_[N];
		for (int k = 0; k <= N / 2; ++k) {
			const int n     = N - 2 * k;
			const T&  sin_n = pows_.sin_[n];
			for (int l = 0; l <= lmax; ++l)
				slab = slab
					   + point_to_ellipse_series::series_coeff<T>(d_h_evo(n, k, l))
						 * pows_.b_[1 - sn + 2 * l] * rho_m * sin_n;
		}
		accum_ = accum_ + slab;
		return slab;
	}

	const T& value() const { return accum_; }

private:
	const EvoBasePowers<T>& pows_;
	T accum_;
};
