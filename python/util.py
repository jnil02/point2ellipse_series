
import sympy as sp
import cache

def rf_half(k: int, n: int):
    """N:th rising factorial of k/2

    :param k: positive integer.
    :param n: positive integer.
    :return:
    """
    r = 1
    s = k
    for i in range(n):
        r = r * s
        s = s + 2
    return sp.Rational(r, 2 ** n)

@cache.ints_cache
def A(n):
    if n==0:
        return 1
    a = sp.S.Zero
    for k in range(n):
        a += sp.binomial(n-1, k) * A(k) * A(n-1-k)
    return a

from math import comb
def E2(n: int) -> int:
    """Euler (secant) number E_{2n}, minimal recursive version."""
    if n == 0:
        return 1
    return sum(((-1) ** (n - j + 1)) * comb(2*n, 2*j) * E2(j) for j in range(n))
