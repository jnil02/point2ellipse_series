#pragma once

/*
 * Point-to-ellipse Fourier and sin-power series expansion coefficients.
 */

namespace point_to_ellipse_series {

// Rational coefficient.
struct rc {
	long num;
	long den;
};

rc c_phi(int n, int k, int l);
rc d_phi(int n, int k, int l);
rc d_phi2(int n, int k, int l);
rc c_sin(int n, int k, int l);
rc d_sin(int n, int k, int l);
rc c_cos(int n, int k, int l);
rc d_cos(int n, int k, int l);
rc c_h(int n, int k, int l);
rc d_h(int n, int k, int l);

}  // namespace point_to_ellipse_series