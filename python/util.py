
""" Various Utility functions for point-to-ellipse series expansions.
"""

# External includes.
from collections.abc import Callable
import sympy as sp
from collections.abc import Hashable


def cantor_pairing_function(k: int, l: int) -> int:
    """Cantor pairing for two non-negative integers."""
    assert k >= 0, "k must be >= 0"
    assert l >= 0, "l must be >= 0"
    return (k + l) * (k + l + 1) // 2 + l

def cantor_tuple(*args: int) -> int:
    """Cantor pairing extended to n non-negative integers."""
    if not args:
        raise ValueError("At least one integer required")
    assert all(isinstance(a, int) and a >= 0 for a in args), \
        "All arguments must be non-negative integers"
    key = args[0]
    for next_val in args[1:]:
        key = cantor_pairing_function(key, next_val)
    return key

def split_args(*args):
    """Split *args into (ints, others) where 'others' are hashable non-ints."""
    ints = []
    others = []
    for a in args:
        if isinstance(a, int):
            ints.append(a)
        elif isinstance(a, Hashable):
            others.append(a)
        else:
            raise TypeError(f"Argument {a!r} is not hashable")
    return ints, tuple(others)

class ints_cache:
    """Caches results of integer-arg + other hashable arg functions."""
    def __init__(self, f):
        self.f = f
        self.cache = {}

    def __call__(self, *args):
        ints, others = split_args(*args)
        if ints:
            int_key = cantor_tuple(*ints)
            key = (int_key,) + others
        else:
            key = others
        if key in self.cache:
            return self.cache[key]
        res = self.f(*args)
        self.cache[key] = res
        return res

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
