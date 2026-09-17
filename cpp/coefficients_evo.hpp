#pragma once

/*
 * Point-to-ellipse inside-evolute series expansion coefficients.
 */

#include "coefficients.hpp"

namespace point_to_ellipse_series {

mpq_class c_phi_evo(int k, int l, int n);
mpq_class d_phi_evo(int k, int l, int n);
mpq_class c_phi_pow_evo(int k, int l, int n, int i);
mpq_class c_cos_phi_evo(int k, int l, int n);
mpq_class d_cos_phi_evo(int k, int l, int n);
mpq_class c_sin_phi_evo(int k, int l, int n);
mpq_class d_sin_phi_evo(int k, int l, int n);
mpq_class c_sin_phi_inv_evo(int k, int l, int n);

mpq_class c_N_evo(int k, int l, int n);
mpq_class cp_evo_nkl(int k, int l, int n);
mpq_class c_h_evo(int k, int l, int n);
mpq_class d_h_evo(int k, int l, int n);

}  // namespace point_to_ellipse_series
