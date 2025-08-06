

"""Example of evaluation of point-to-ellipse Fourier series.

The earth ellipse is used as an example.
All non-integer evaluations are done with multi-precision arithmetics.

An arbitrary point in latitude and altitude is chosen for reference and
converted to an ECEF (Cartesian) coordinate. Subsequently, all the series with
coefficients are computed up to a specified order (truncated) and evaluated.
Finally, their values and the original reference values are printed for
comparison.
"""

# Multiprecision arithmetics.
from mpmath import mp

# Set the number of digits of precision for the calculations.
mp.dps = 50

# Utility stuff.
from ellipse import mp_ellipse_to_cartesian, mp_e2, mp_a
from symbols import varrho, psi, sin_psi, cos_psi, e2

# The actual series.
import series

if __name__ == "__main__":
    # Number of terms in each series.
    max_order = 7

    # Multi-precision reference/input values.
    mp_phi = mp.mpf(43.099999999999504)/mp.mpf(180.) * mp.pi  # "Latitude".
    mp_h = mp.mpf(10000.)  # "Altitude".
    mp_x, mp_y = mp_ellipse_to_cartesian(mp_phi, mp_h)  # Cartesian coordinate.
    mp_psi = mp.atan2(mp_y , mp_x)  # Polar angle.
    mp_rho = mp.sqrt(mp_x * mp_x + mp_y * mp_y)  # Radius.
    mp_varrho = mp_a / mp.sqrt(mp_x * mp_x + mp_y * mp_y)

    # For each series quentity:
    # 1. compute and print reference value,
    # 2. compute symbolic series up to given order,
    # 3. and substitute numerica values in series and print numeric value.

    print("phi - psi (ref/sin-mul/sin-pow):")
    print(mp_phi-mp_psi)
    expr = series.phi_in_sin_mul(max_order, max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(psi, mp_psi)
          .subs(varrho, mp_varrho))
    expr = series.phi_in_sin_pow(max_order, max_order) * sin_psi * cos_psi
    print(expr.subs(e2, mp_e2)
          .subs(sin_psi, mp.sin(mp_psi))
          .subs(cos_psi, mp.cos(mp_psi))
          .subs(varrho, mp_varrho))

    print("psi in sin pow with cos-factor integrated")
    expr = series.phi_in_sin_pow2(max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(sin_psi, mp.sin(mp_psi))
          .subs(cos_psi, mp.cos(mp_psi))
          .subs(varrho, mp_varrho))

    print("sin(phi)/sin(psi) (ref/cos-mul/sin-pow):")
    if (mp_psi != 0):
        print(mp.sin(mp_phi) / mp.sin(mp_psi))
    else:
        print("NaN")
    expr = series.sin_phi_in_sin_pow(max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(sin_psi, mp.sin(mp_psi))
          .subs(varrho, mp_varrho) + 1)
    expr = series.sin_phi_in_cos_mul(max_order, max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(psi, mp_psi)
          .subs(varrho, mp_varrho) + 1)
    expr = series.sin_phi_in_sin_pow2(max_order, max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(sin_psi, mp.sin(mp_psi))
          .subs(varrho, mp_varrho))

    print("cos(phi)/cos(psi) (ref/cos-mul/sin-pow):")
    print(mp.cos(mp_phi) / mp.cos(mp_psi))
    expr = series.cos_phi_in_sin_pow(max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(sin_psi, mp.sin(mp_psi))
          .subs(cos_psi, mp.cos(mp_psi))
          .subs(varrho, mp_varrho) + 1)
    expr = series.cos_phi_in_cos_mul(max_order, max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(psi, mp_psi)
          .subs(varrho, mp_varrho) + 1)
    expr = series.cos_phi_in_sin_pow2(max_order, max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(sin_psi, mp.sin(mp_psi))
          .subs(varrho, mp_varrho))

    print("(h+a-rho)/a (ref/cos-mul/sin-pow):")
    print((mp_h + mp_a - mp.sqrt(mp_x * mp_x + mp_y * mp_y)) / mp_a)
    expr = series.h_in_sin_pow(max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(sin_psi, mp.sin(mp_psi))
          .subs(varrho, mp_varrho))
    expr = series.h_in_cos_mul(max_order, max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(psi, mp_psi)
          .subs(varrho, mp_varrho))

    print("h:  (ref/sin-pow)")
    print(mp_h)
    expr = series.h_in_sin_pow(max_order, max_order)
    print(expr.subs(e2, mp_e2)
          .subs(sin_psi, mp.sin(mp_psi))
          .subs(varrho, mp_varrho) * mp_a + mp_rho - mp_a)

