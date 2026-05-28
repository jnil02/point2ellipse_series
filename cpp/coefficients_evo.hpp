#pragma once

/*
 * Point-to-ellipse inside-evolute series expansion coefficients.
 */

#include "coefficients.hpp"

namespace point_to_ellipse_series {

mpq_class d_phi_evo(int n, int k, int l);
mpq_class d_sin_phi_evo(int n, int k, int l);
mpq_class d_cos_phi_evo(int n, int k, int l);
mpq_class c_phi_evo(int n, int k, int l);
mpq_class c_phi_pow_evo(int n, int k, int l, int i);
mpq_class c_sin_phi_evo(int n, int k, int l);
mpq_class c_cos_phi_evo(int n, int k, int l);
mpq_class c_sin_phi_inv_evo(int n, int k, int l);

// TODO(JO) Temporary. These are intermediate coefficients and should be
//  removed once the final coefficients are in place.
mpq_class a_mr(int m, int r);
mpq_class B_rt(int r, int t);
mpq_class C_mt(int m, int t);
mpq_class R(int n, int k, int l, int i);
mpq_class B_p(int n, int k, int p);
mpq_class cp_evo_nkl(int n, int k, int l);
mpq_class c_h_evo(int n, int k, int l);
mpq_class d_h_evo(int n, int k, int l);

// FIXME(JO) Exposed for tests. Put in different header.
mpq_class c_phi_pow_evo_se(int n, int k, int l, int i);
mpq_class c_phi_pow_evo_se2(int n, int k, int l, int i);

}  // namespace point_to_ellipse_series
