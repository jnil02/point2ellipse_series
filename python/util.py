
""" Various Utility functions for point-to-ellipse series expansions.
"""

# External includes.
from collections.abc import Callable
import sympy as sp


def sin_pow_to_cos_mul(n: int, k: int, l: int, d_nkl: Callable[[int, int, int], sp.core.Expr],
                       n_min: int, k_pp: int) -> sp.core.Rational | sp.core.Integer:
    """Fourier multiple-angle cos series coefficient from sin-power series.

    :param n: sin-multiple.
    :param k: rho power.
    :param l: e³ power.
    :param d_nkl: Sin-power series coefficients.
    :param n_min: Lowest sin-power.
    :param k_pp: Maximum e² power offset in sin-power series.
    :return: Fourier multiple-angle cos series coefficient of
    """
    c_nkl = sp.Integer(0)
    for i in range(max(max(n, n_min), l - k - k_pp), l + 1):
        tmp = d_nkl(i, k, l) * sp.Rational(2, 2 ** (2 * i)) * (-1) ** n * sp.binomial(2 * i, i - n)
        c_nkl = c_nkl + (tmp if n != 0 else tmp / 2)
    return c_nkl
