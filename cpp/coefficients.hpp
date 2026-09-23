#pragma once

/*
 * Point-to-ellipse Fourier and sin-power series expansion coefficients (far-field).
 */

#include <gmpxx.h>

namespace point_to_ellipse_series {

/// Fourier-series coefficient for φ−ψ.
mpq_class c_phi(int n, int k, int l);
/// sin-power-series coefficient for φ−ψ.
mpq_class d_phi(int n, int k, int l);
/// sin-power variant of d_phi with the cos·sin factor folded in.
mpq_class d_phi2(int n, int k, int l);
/// Fourier-series coefficient for sin(φ)/sin(ψ)−1.
mpq_class c_sin(int n, int k, int l);
/// sin-power-series coefficient for sin(φ)/sin(ψ)−1.
mpq_class d_sin(int n, int k, int l);
/// Fourier-series coefficient for cos(φ)/cos(ψ)−1.
mpq_class c_cos(int n, int k, int l);
/// sin-power-series coefficient for cos(φ)/cos(ψ)−1.
mpq_class d_cos(int n, int k, int l);
/// Fourier-series coefficient for the height term (h+a−ρ)/a.
mpq_class c_h(int n, int k, int l);
/// sin-power-series coefficient for the height term (h+a−ρ)/a.
mpq_class d_h(int n, int k, int l);

}  // namespace point_to_ellipse_series
