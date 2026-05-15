#pragma once

/*
 * Point-to-ellipse inside-evolute series expansion coefficients.
 */

#include "coefficients.hpp"

namespace point_to_ellipse_series {

rc d_phi_evo(int n, int k, int l);
rc c_phi_evo(int n, int k, int l);
rc c_phi_pow_evo(int n, int k, int l, int i);
rc c_sin_phi_evo(int n, int k, int l);
rc c_cos_phi_evo(int n, int k, int l);
rc c_sin_phi_inv_evo(int n, int k, int l);

// TODO(JO) Temporary. These are intermediate coefficients and should be
//  removed once the final coefficients are in place.
rc a_mr(int m, int r);
rc B_rt(int r, int t);
rc C_mt(int m, int t);
rc R(int n, int k, int l, int i);
rc B_p(int n, int k, int p);
rc cp_evo_nkl(int n, int k, int l);
rc ch_evo(int n, int k, int l);
rc dh_evo(int n, int k, int l);

}  // namespace point_to_ellipse_series
