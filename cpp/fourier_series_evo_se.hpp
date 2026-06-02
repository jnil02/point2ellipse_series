#pragma once

#include "fourier_series_evo.hpp"
#include "fourier_series_se.hpp"

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_evo_sin_pow_dense(int N, int K) {
	return phi_evo_sin_pow_dense<Expression>(N, K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_evo_sin_pow_dense_m(int K) {
	return phi_evo_sin_pow_dense_m<Expression>(K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression phi_evo_sin_pow(int N, int K) {
	return phi_evo_sin_pow<Expression>(N, K, sin_psi, cos_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression sin_phi_evo_dense(int N, int K) {
	return sin_phi_evo_dense<Expression>(N, K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression sin_phi_evo_dense_m(int K) {
	return sin_phi_evo_dense_m<Expression>(K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression cos_phi_evo_dense(int N, int K) {
	return cos_phi_evo_dense<Expression>(N, K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression cos_phi_evo_dense_m(int N) {
	return cos_phi_evo_dense_m<Expression>(N, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_a_evo(int N, int K) {
	return h_a_evo<Expression>(N, K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_a_evo_dense(int N, int K) {
	return h_a_evo_dense<Expression>(N, K, sin_psi, rho_ae2, b_a);
}

/** Symbolic convenience overload: returns Expression using the global symbolic variables. */
inline Expression h_a_evo_dense_m(int K) {
	return h_a_evo_dense_m<Expression>(K, sin_psi, rho_ae2, b_a);
}
