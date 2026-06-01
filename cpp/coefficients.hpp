#pragma once

/*
 * Point-to-ellipse Fourier and sin-power series expansion coefficients.
 */

#include <gmpxx.h>

namespace point_to_ellipse_series {

mpq_class c_phi(int n, int k, int l);
mpq_class d_phi(int n, int k, int l);
mpq_class d_phi2(int n, int k, int l);
mpq_class c_sin(int n, int k, int l);
mpq_class d_sin(int n, int k, int l);
mpq_class c_cos(int n, int k, int l);
mpq_class d_cos(int n, int k, int l);
mpq_class c_h(int n, int k, int l);
mpq_class d_h(int n, int k, int l);

}  // namespace point_to_ellipse_series
