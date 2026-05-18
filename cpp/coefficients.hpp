#pragma once

/*
 * Point-to-ellipse Fourier and sin-power series expansion coefficients.
 */

namespace point_to_ellipse_series {

// Rational coefficient.
// __int128 is used instead of long to avoid overflow at higher truncation
// orders. long, with its 9 decimal digits, overflows around order 13 while
// __int128 extends the range to ~39 decimal digits.
struct rc {
	__int128 num;
	__int128 den;
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
