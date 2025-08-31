#pragma once

#include <utility>  // std::pair
#include <mpreal.h>

using mpfr::mpreal;
using mpfr::abs;
using mpfr::sin;
using mpfr::cos;
using mpfr::atan;
using mpfr::sqrt;
using mpfr::cbrt;
using mpfr::const_pi;

// Set default precision once for your app (like mp.dps but in bits).
inline void set_precision_bits(mpfr_prec_t bits) {
	mpreal::set_default_prec(bits);
}

// WGS84 constants. Lazy evaluation such that set precision is respected.
// The literal values are defined by the standard.
inline const mpreal& mp_a()  { static const mpreal v = mpreal("6378137.0");                               return v; } // Semi-major axis.
inline const mpreal& mp_f()  { static const mpreal v = mpreal(1) / mpreal("298.257223563");               return v; } // Flattening.
inline const mpreal& mp_b()  { static const mpreal v = mp_a() - mp_f() * mp_a();                          return v; } // Semi-minor axis.
inline const mpreal& mp_e2() { static const mpreal v = mpreal(1) - (mp_b() * mp_b()) / (mp_a() * mp_a()); return v; } // Eccentricity squared.

// Geodetic (lat,alt) to 2D Cartesian (x,y).
inline std::pair<mpreal, mpreal>
mp_ellipse_to_cartesian(const mpreal& phi, const mpreal& h) {
	const mpreal s = sin(phi);
	const mpreal c = cos(phi);
	const mpreal N = mp_a() / sqrt(mpreal(1) - mp_e2() * s * s);
	const mpreal x = (N + h) * c;
	const mpreal y = ((mpreal(1) - mp_e2()) * N + h) * s;
	return {x, y};
}

// Cartesian to elliptical coordinate transformation by Vermeille (2011).
inline std::pair<mpreal, mpreal>
mp_cartesian_to_ellipse(const mpreal& x, const mpreal& z) {
	// Derived constants. Cached on first call with the precision at that time.
	static const mpreal a2          = mp_a() * mp_a();
	static const mpreal b2          = mp_b() * mp_b();
	static const mpreal e4          = mp_e2() * mp_e2();
	static const mpreal b2oa2       = b2 / a2;
	static const mpreal a3          = a2 * mp_a();
	static const mpreal a4          = a2 * a2;
	static const mpreal inv_a2      = mpreal(1) / a2;
	static const mpreal b_over_a2   = mp_b() / a2;
	static const mpreal b2_over_a4  = b2 / a4;
	static const mpreal inv_6       = mpreal(1) / mpreal(6);
	static const mpreal e2_over_2   = mp_e2() / mpreal(2);
	static const mpreal e2b_over_a3 = mp_e2() * mp_b() / a3;
	static const mpreal PI_over_6   = const_pi() / mpreal(6);
	static const mpreal two_thirds  = mpreal(2) / mpreal(3);
	static const mpreal half        = mpreal(0.5);
	static const mpreal b_over_a    = mp_b() / mp_a();
	static const mpreal b_over_e    = mp_b() / sqrt(mp_e2());
	static const mpreal e2_over_b   = mp_e2() / mp_b();

	const mpreal t2  = x * x;
	const mpreal t   = sqrt(t2);
	const mpreal p   = t2 * inv_a2;
	const mpreal q   = z * z * b2_over_a4;
	const mpreal r   = (p + q - e4) * inv_6;
	const mpreal r38 = r * r * r * mpreal(8);
	const mpreal ev  = r38 + e4 * p * q;

	if (ev > 0) {
		const mpreal s = e2b_over_a3 * abs(z) * t;
		const mpreal sqrt_ev = sqrt(ev);
		const mpreal c1 = sqrt_ev + s;
		const mpreal c2 = sqrt_ev - s;

		const mpreal u = r + half * cbrt(c1 * c1) + half * cbrt(c2 * c2);
		const mpreal v = sqrt(u * u + e4 * q);
		const mpreal w = e2_over_2 * (u + v - q) / v;
		const mpreal k = (u + v) / (sqrt(w * w + u + v) + w);
		const mpreal D = k / (k + mp_e2()) * t;
		const mpreal d = sqrt(D * D + z * z);

		const mpreal h   = (k - b2oa2) * d / k;
		const mpreal lat = mpreal(2) * atan(z / (d + D));
		return {lat, h};
	} else if (q != 0) {
		const mpreal s  = e2b_over_a3 * abs(z) * t;
		const mpreal up = two_thirds * atan(s / (sqrt(-ev) + sqrt(-r38)));
		const mpreal u  = mpreal(-4) * r * sin(up) * cos(PI_over_6 + up);

		const mpreal v = sqrt(u * u + e4 * q);
		const mpreal w = e2_over_2 * (u + v - q) / v;
		const mpreal k = (u + v) / (sqrt(w * w + u + v) + w);
		const mpreal D = k * t / (k + mp_e2());
		const mpreal d = sqrt(D * D + z * z);

		const mpreal h   = (k - b2oa2) * d / k;
		const mpreal lat = mpreal(2) * atan(z / (d + D));
		return {lat, h};
	} else {
		const mpreal h   = -b_over_e * sqrt(mp_e2() - p);
		const mpreal lat = mpreal(2) * atan(sqrt(e4 - p) / (-e2_over_b * h + b_over_a * sqrt(p)) );
		return {lat, h};
	}
}
