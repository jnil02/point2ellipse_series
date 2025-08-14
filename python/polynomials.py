
""" Polynomials for series of power of series handling.
"""

# External includes.
import sympy as sp

# Internal includes.
import util

@util.iis_cache
def partial_ordinary_bell_polynomial(k : int, i : int, a : str) -> sp.core.Expr:
    """Partial ordinary Bell polynomial.

    :param k: Power index of polynomial.
    :param i: Sum index of polynomial.
    :param a: Base name for the variables of which the polynomial is a function of.
    :return: Sympy expression for partial ordinary Bell polynomial.
    """
    # Recursion return point.
    if i == 0:
        return sp.Integer(1) if k == 0 else sp.Integer(0)
    # Recursive computation of partial ordinary Bell polynomial.
    tmp = sp.Integer(0)
    for j in range(1, k - i + 2):
        tmp = tmp + sp.symbols(a + "_" + str(j)) * partial_ordinary_bell_polynomial(k - j, i - 1, a)
    return sp.expand(tmp)

@util.iis_cache
def power_of_power_series_coefficient_polynomial(n : int, i : int, a : str) -> sp.core.Expr:
    """Power of power-series series coefficient polynomial.

    Polynomial for n:th series coefficient of i:th power of an infinite series.

    Let p be a power series

    p = \sum_{n=0}^\infty a_n x^n

    where a_0 not equal 0. Then

    p^j = \sum_{n=0}^\infty b_{n,j} x^n

    where

    b_{0,j} = a_0^j
    b_{n,j} = 1/(n b_{0,j}) \sum_{k=1}^n (kj - n + k) a_k b_{n-k,j}

    The recurrence relation for b_{n,j} defines a polynomial

    b_{n,i} = P_{n,i}(a_0,...,a_n)

    :param n: The number of the coefficient to return the polynomial for.
    :param i: The power of the original series.
    :param a: The base symbol name of the original coefficients.
    :return: The n:th coefficient of the i:th power of the power series.
    """
    x0 = sp.symbols(a +"_0")
    # Recursion return point.
    if n==0:
        return x0 ** i
    # Recursive evaluation of polynomial.
    clj = 0
    for k in range(1, n+1):
        clj = (clj + (k * i - n + k) * sp.symbols(a + "_" + str(k)) *
               power_of_power_series_coefficient_polynomial(n - k, i, a))
    return clj / (n * x0)
