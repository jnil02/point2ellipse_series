#pragma once

#include "fourier_series_evo.hpp"
#include "fourier_series_se.hpp"

using point_to_ellipse_series::phi_evo_sparse;
using point_to_ellipse_series::phi_evo_dense;
using point_to_ellipse_series::sin_phi_evo_sparse;
using point_to_ellipse_series::sin_phi_evo_dense;
using point_to_ellipse_series::cos_phi_evo_sparse;
using point_to_ellipse_series::cos_phi_evo_dense;
using point_to_ellipse_series::h_evo_sparse;
using point_to_ellipse_series::h_evo_dense;

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
