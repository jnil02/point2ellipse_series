#pragma once

#include "fourier_series_evo.hpp"
#include "fourier_series_se.hpp"

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_evo_sparse(int L, int K) {
	return phi_evo_sparse<Expression>(L, K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_evo_dense(int K) {
	return phi_evo_dense<Expression>(K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression sin_phi_evo_dense(int K) {
	return sin_phi_evo_dense<Expression>(K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression cos_phi_evo_dense(int K) {
	return cos_phi_evo_dense<Expression>(K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_evo_sparse(int L, int K) {
	return h_evo_sparse<Expression>(L, K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_evo_dense(int K) {
	return h_evo_dense<Expression>(K, sin_psi, rho_ae2, b_a);
}
