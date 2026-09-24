#pragma once

/*
 * Point-to-ellipse inside-evolute series expansion coefficients.
 */

#include <gmpxx.h>

namespace point_to_ellipse_series {

/// Sparse φ−π/2 coefficient.
mpq_class c_phi_evo(int k, int l, int n);
/// Dense φ−π/2 coefficient.
mpq_class d_phi_evo(int k, int l, int n);
/// Sparse (φ−π/2)^i coefficient.
mpq_class c_phi_pow_evo(int k, int l, int n, int i);
/// Sparse cos(φ)/cos(ψ) coefficient.
mpq_class c_cos_phi_evo(int k, int l, int n);
/// Dense cos(φ)/cos(ψ) coefficient.
mpq_class d_cos_phi_evo(int k, int l, int n);
/// Sparse sin(φ) coefficient.
mpq_class c_sin_phi_evo(int k, int l, int n);
/// Dense sin(φ) coefficient.
mpq_class d_sin_phi_evo(int k, int l, int n);
/// Sparse sin(ψ)/sin(φ) coefficient.
mpq_class c_sin_phi_inv_evo(int k, int l, int n);

/// Sparse ε²·N/a coefficient (N the normal radius of curvature).
mpq_class c_N_evo(int k, int l, int n);
/// Sparse sin(ψ)/sin(φ) contribution to h (ρ·sin(ψ)/a pulled out).
mpq_class cp_evo(int k, int l, int n);
/// Sparse h/a coefficient (ρ·sin(ψ)/a pulled out).
mpq_class c_h_evo(int k, int l, int n);
/// Dense h/a coefficient (ρ·sin(ψ)/a pulled out).
mpq_class d_h_evo(int k, int l, int n);

}  // namespace point_to_ellipse_series