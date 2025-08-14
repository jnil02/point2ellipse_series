
""" Various Utility functions for point-to-ellipse series expansions.
"""

# External includes.
from collections.abc import Callable
import sympy as sp

def cantor_pairing_function(k : int, l : int) -> int:
    """Cantor pairing function.

    :param k: A positive integer.
    :param l: Another positive integer.
    :return: A unique bijection from k and l to another positive integer.
    """
    return int( (k + l) * (k + l + 1) / 2 + l )

class ii_cache:
    """Dual integer argument function cache decorator.
    """
    def __init__(self, f):
        self.f = f
        self.cache = {}
    def __call__(self, i : int, j : int):
        key = cantor_pairing_function(i, j)
        if key in self.cache:
            return self.cache[key]
        res = self.f(i, j)
        self.cache[key] = res
        return res

class iii_cache:
    """Triple integer argument function cache decorator.
    """
    def __init__(self, f):
        self.f = f
        self.cache = {}
    def __call__(self, i : int, j : int, k : int):
        key = cantor_pairing_function(cantor_pairing_function(i, j), k)
        if key in self.cache:
            return self.cache[key]
        res = self.f(i, j, k)
        self.cache[key] = res
        return res

class iiii_cache:
    """Quad integer argument function cache decorator.
    """
    def __init__(self, f):
        self.f = f
        self.cache = {}
    def __call__(self, i : int, j : int, k : int, l : int):
        key = cantor_pairing_function(cantor_pairing_function(i, j), cantor_pairing_function(k, l))
        if key in self.cache:
            return self.cache[key]
        res = self.f(i, j, k, l)
        self.cache[key] = res
        return res

class iis_cache:
    """Dual integer and string argument function cache decorator.
    """
    def __init__(self, f):
        self.f = f
        self.cache = {}
    def __call__(self, i : int, j : int, s : str):
        if s not in self.cache:
            self.cache[s] = {}
        scache = self.cache[s]
        key = cantor_pairing_function(i, j)
        if key in scache:
            return scache[key]
        res = self.f(i, j, s)
        scache[key] = res
        return res


def sin_pow_to_cos_mul(n : int, k : int, l : int, d_nkl : Callable[[int, int, int], sp.core.Expr],
                       n_min : int, k_pp : int) -> sp.core.Rational | sp.core.Integer:
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
