
import sympy as sp


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
