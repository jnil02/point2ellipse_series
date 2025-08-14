
"""Series expansions for point-to-ellipse relations.

For further details see
Nilsson, John-Olof. Point-to-ellipse Fourier series. DOI:
https://doi.org/10.48550/arXiv.2506.XXXXX.
"""

import sympy as sp

from symbols import varrho, psi, sin_psi, cos_psi, e2
from coefficients import c_phi, d_phi, d_phi2, c_sin, c_cos, d_phi_pow, d_cos, d_sin, c_h, d_h

def phi_in_sin_pow(n_max: int, k_max: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of phi-psi up to given order.

    :param n_max: Maximum sin²-power.
    :param k_max: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, n_max+1):
        for k in range(1, k_max+1):
            for l in range(max(k,n+1), k+n + 1):
                d = d + d_phi(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d

def phi_in_sin_pow2(n_max: int, k_max: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of phi-psi up to given order.

    Series with the proceeding cos incorporated into the series.
    This gives worse convergence.

    :param n_max: Maximum sin²-power.
    :param k_max: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, n_max+1):
        for k in range(1, k_max+1):
            for l in range(k, k+n + 1):
                d = d + d_phi2(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n + 1)
    return d

def phi_in_sin_mul(n_max: int, k_max: int, l_max: int) -> sp.core.Expr:
    """Symbolic sin (Fourier) series expansion of phi-psi up to given order.

    :param n_max: Maximum sin-multiple.
    :param k_max: Maximum varrho power.
    :param l_max: Maxiumu e² power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(1, n_max+1):
        for k in range(1, k_max+1):
            for l in range(max(n,k), l_max+1):
                d = d + c_phi(n,k,l) * e2 ** l * varrho ** k * sp.sin(2 * n * psi)
    return d

def phi_pow_in_sin_pow(i: int, n_max: int, k_max: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of (phi-psi)^i up to given order.

    :param i: The power exponent.
    :param n_max: Maximum sin²-power.
    :param k_max: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, n_max+1):
        for k in range(i, k_max+1):
            for l in range(max(k,n+i), n+k+1):
                d = d + d_phi_pow(n, k, l, i) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d * cos_psi ** i * sin_psi ** i

def h_in_cos_mul(n_max: int, k_max: int, l_max: int) -> sp.core.Expr:
    """Symbolic cos (Fourier) series expansion of (h+a-rho)/a up to given order.

    :param n_max: Maximum cos-multiple.
    :param k_max: Maximum varrho power.
    :param l_max: Maxiumu e² power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, n_max+1):
        for k in range(k_max+1):
            for l in range(max(n,k+1), l_max+1):
                d = d + c_h(n,k,l) * e2 ** l * varrho ** k * sp.cos(2 * n * psi)
    return d

def h_in_sin_pow(n_max: int, k_max: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of (h+a-rho)/a up to given order.

    :param n_max: Maximum sin²-power.
    :param k_max: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(1, n_max+1):
        for k in range(0, k_max+1):
            for l in range(max(k+1,n), n+k+1):
                d = d + d_h(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d

def sin_phi_in_cos_mul(n_max: int, k_max: int, l_max: int) -> sp.core.Expr:
    """Symbolic cos (Fourier) series expansion of sin(phi)/sin(psi) - 1 up to given order.

    :param n_max: Maximum cos-multiple.
    :param k_max: Maximum varrho power.
    :param l_max: Maxiumu e² power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, n_max+1):
        for k in range(1, k_max+1):
            for l in range(max(n,k), l_max+1):
                d = d + c_sin(n,k,l) * e2 ** l * varrho ** k * sp.cos(2 * n * psi)
    return d

def sin_phi_in_sin_pow(n_max: int, k_max: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of sin(phi)/sin(psi) - 1 up to given order.

    :param n_max: Maximum sin²-power.
    :param k_max: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, n_max+1):
        for k in range(1, k_max+1):
            for l in range(max(k,n), n+k+1):
                d = d + d_sin(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d

def cos_phi_in_cos_mul(n_max: int, k_max: int, l_max: int) -> sp.core.Expr:
    """Symbolic cos (Fourier) series expansion of cos(phi)/cos(psi) - 1 up to given order.

    :param n_max: Maximum cos-multiple.
    :param k_max: Maximum varrho power.
    :param l_max: Maxiumu e² power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, n_max+1):
        for k in range(1, k_max+1):
            for l in range(max(n,k), l_max+1):
                d = d + c_cos(n,k,l) * e2 ** l * varrho ** k * sp.cos(2 * n * psi)
    return d

def cos_phi_in_sin_pow(n_max: int, k_max: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of cos(phi)/cos(psi) - 1 up to given order.

    :param n_max: Maximum sin²-power.
    :param k_max: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, n_max+1):
        for k in range(1, k_max+1):
            for l in range(max(k,n), n+k):
                d = d + d_cos(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d

def sigma(j_max: int, delta: sp.core.Expr) -> sp.core.Expr:
    """ Symbolic "sigma" series up to given order.

    :param j_max: Maximum delta power.
    :param delta: Expression for delta.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for j in range(0, j_max, 2):
        j2 = j / 2
        d = d + sp.Rational(sp.Integer(-1) ** j2, sp.factorial(j)) * delta ** j2
    return d

def tau(j_max: int, omega: sp.core.Expr, delta: sp.core.Expr):
    """ Symbolic "tau" series up to given order.

    :param j_max: Maximum delta power.
    :param omega: Expression for "omega".
    :param delta: Expression for delta.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for j in range(1, j_max, 2):
        j2 = (j - 1) / 2
        d = d + sp.Rational(sp.Integer(-1) ** j2, sp.factorial(j)) * delta ** j2
    return omega * d

def sin_phi_in_sin_pow2(n_max: int, k_max: int, j_max: int):
    """Symbolic sin-power series expansion of sin(phi)/sin(psi) up to given order.

    Using sigma and tau series.

    :param n_max: Maximum sin²-power.
    :param k_max: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    sin_psi2 = sin_psi ** 2
    o = phi_in_sin_pow(n_max, k_max)
    d = sin_psi2 * (sp.Integer(1) - sin_psi2) * o ** 2
    s = sigma(j_max, d)
    t = tau(j_max, o, d)
    return s + (sp.Integer(1) - sin_psi2) * t

def cos_phi_in_sin_pow2(n_max: int, k_max: int, j_max: int):
    """Symbolic sin-power series expansion of sin(phi)/sin(psi) up to given order.

    Using sigma and tau series.

    :param n_max: Maximum sin²-power.
    :param k_max: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    sin_psi2 = sin_psi ** 2
    o = phi_in_sin_pow(n_max, k_max)
    d = sin_psi2 * (sp.Integer(1) - sin_psi2) * o ** 2
    s = sigma(j_max, d)
    t = tau(j_max, o, d)
    return s - sin_psi2 * t
