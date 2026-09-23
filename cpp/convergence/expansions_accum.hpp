#pragma once

/*
 * Incremental (slab-by-slab) accumulators for the inside-evolute series.
 *
 * Evaluating phi_evo_dense or h_evo_dense at each truncation order N in a loop
 * re-sums all prior slabs, giving O(N^2) work per (psi, rho) point and O(N^3)
 * total.  The classes here decouple slab addition so each new order costs only
 * the m=N contribution — O(N^2) total.  Each addOrder(N) reproduces exactly the
 * k=N term of the corresponding dense series (same coefficient, same powers).
 *
 * Usage:
 *
 *   EvoBasePowers<mpreal> pows(sin_psi, rho_ae2, b_a, MAX_ORDER);
 *   PhiEvoAccum<mpreal>   phi_acc(pows);
 *   HEvoAccum<mpreal>    h_acc(pows);
 *   for (int N = 1; N <= MAX_ORDER; ++N) {
 *       mpreal cm = phi_acc.addOrder(N);  // returns slab C_m; also updates value()
 *       h_acc.addOrder(N);                // h_evo_dense starts at k=2; N<2 adds 0
 *       use(phi_acc.value(), h_acc.value());
 *   }
 *
 * EvoBasePowers holds power tables shared across all accumulators for the same
 * (psi, rho) point; adding new series requires only a new accumulator class.
 */

#include <vector>

#include "coefficients_evo.hpp"
#include "series_traits.hpp"

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
// PhiEvoAccum — incremental accumulator for phi_evo_dense
// ---------------------------------------------------------------------------

template<typename T>
class PhiEvoAccum {
public:
	explicit PhiEvoAccum(const EvoBasePowers<T>& pows)
			: pows_(pows), accum_(T(0)) {}

	// Add the m=N slab (mirrors the k=N term of phi_evo_dense).
	// Returns the slab value so callers can inspect C_m without a separate class.
	T addOrder(int N) {
		T slab(0);
		const int s = (N + 1) % 2;
		const int r = (N - 1) / 2;
		for (int l = 0; l <= r; ++l)
			for (int n = 0; n <= r; ++n)
				slab = slab
					   + point_to_ellipse_series::to<T>(
						point_to_ellipse_series::d_phi_evo(N, l, n))
						 * pows_.sin_[2 * l] * pows_.b_[2 * n];
		slab = slab * pows_.rho_[N] * pows_.sin_[s] * pows_.b_[s + 1];
		accum_ = accum_ + slab;
		return slab;
	}

	const T& value() const { return accum_; }

private:
	const EvoBasePowers<T>& pows_;
	T accum_;
};

// ---------------------------------------------------------------------------
// HEvoAccum — incremental accumulator for h_evo_dense
// ---------------------------------------------------------------------------

template<typename T>
class HEvoAccum {
public:
	explicit HEvoAccum(const EvoBasePowers<T>& pows)
			: pows_(pows), accum_(T(0)) {}

	// Add the m=N slab (mirrors the k=N term of h_evo_dense); zero for N<2.
	// Returns the slab value.
	T addOrder(int N) {
		T slab(0);
		if (N >= 2) {
			const int p = N % 2;
			for (int l = 0; l <= N / 2; ++l)
				for (int n = 0; n <= N / 2; ++n)
					slab = slab
						   + point_to_ellipse_series::to<T>(
							point_to_ellipse_series::d_h_evo(N, l, n))
							 * pows_.b_[2 * n] * pows_.sin_[2 * l];
			slab = slab * pows_.b_[1 + p] * pows_.sin_[p] * pows_.rho_[N];
		}
		accum_ = accum_ + slab;
		return slab;
	}

	const T& value() const { return accum_; }

private:
	const EvoBasePowers<T>& pows_;
	T accum_;
};
