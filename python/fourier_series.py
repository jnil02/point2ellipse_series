
"""
Symbolic series expansions for the point-to-ellipse relation.
"""

import sympy as sp

from symbols import varrho, rho_ae2, psi, sin_psi, cos_psi, e2, b_a
from coefficients import c_phi, d_phi, d_phi2, c_sin, c_cos, d_phi_pow, d_cos, d_sin, c_h, d_h, d_phi_evo_dense, c_phi_evo, d_phi_pow_evo

def phi_in_sin_pow(N: int, K: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of phi-psi up to given order.

    :param N: Maximum sin²-power.
    :param K: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, N+1):
        for k in range(1, K+1):
            for l in range(max(k,n+1), k+n + 1):
                d = d + d_phi(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d

def phi_in_sin_pow2(N: int, K: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of phi-psi up to given order.

    Series with the proceeding cos incorporated into the series.
    This gives worse convergence.

    :param N: Maximum sin²-power.
    :param K: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, N+1):
        for k in range(1, K+1):
            for l in range(k, k+n + 1):
                d = d + d_phi2(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n + 1)
    return d

def phi_in_sin_mul(N: int, K: int, L: int) -> sp.core.Expr:
    """Symbolic sin (Fourier) series expansion of phi-psi up to given order.

    :param N: Maximum sin-multiple.
    :param K: Maximum varrho power.
    :param L: Maxiumu e² power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(1, N+1):
        for k in range(1, K+1):
            for l in range(max(n,k), L+1):
                d = d + c_phi(n,k,l) * e2 ** l * varrho ** k * sp.sin(2 * n * psi)
    return d

def phi_pow_in_sin_pow(i: int, N: int, K: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of (phi-psi)^i up to given order.

    :param i: The power exponent.
    :param N: Maximum sin²-power.
    :param K: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, N+1):
        for k in range(i, K+1):
            for l in range(max(k,n+i), n+k+1):
                d = d + d_phi_pow(n, k, l, i) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d * cos_psi ** i * sin_psi ** i

def h_in_cos_mul(N: int, K: int, L: int) -> sp.core.Expr:
    """Symbolic cos (Fourier) series expansion of (h+a-rho)/a up to given order.

    :param N: Maximum cos-multiple.
    :param K: Maximum varrho power.
    :param L: Maxiumu e² power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, N+1):
        for k in range(K+1):
            for l in range(max(n,k+1), L+1):
                d = d + c_h(n,k,l) * e2 ** l * varrho ** k * sp.cos(2 * n * psi)
    return d

def h_in_sin_pow(N: int, K: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of (h+a-rho)/a up to given order.

    :param N: Maximum sin²-power.
    :param K: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(1, N+1):
        for k in range(0, K+1):
            for l in range(max(k+1,n), n+k+1):
                d = d + d_h(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d

def sin_phi_in_cos_mul(N: int, K: int, L: int) -> sp.core.Expr:
    """Symbolic cos (Fourier) series expansion of sin(phi)/sin(psi) - 1 up to given order.

    :param N: Maximum cos-multiple.
    :param K: Maximum varrho power.
    :param L: Maxiumu e² power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, N+1):
        for k in range(1, K+1):
            for l in range(max(n,k), L+1):
                d = d + c_sin(n,k,l) * e2 ** l * varrho ** k * sp.cos(2 * n * psi)
    return d

def sin_phi_in_sin_pow(N: int, K: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of sin(phi)/sin(psi) - 1 up to given order.

    :param N: Maximum sin²-power.
    :param K: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, N+1):
        for k in range(1, K+1):
            for l in range(max(k,n), n+k+1):
                d = d + d_sin(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d

def cos_phi_in_cos_mul(N: int, K: int, L: int) -> sp.core.Expr:
    """Symbolic cos (Fourier) series expansion of cos(phi)/cos(psi) - 1 up to given order.

    :param N: Maximum cos-multiple.
    :param K: Maximum varrho power.
    :param L: Maxiumu e² power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, N+1):
        for k in range(1, K+1):
            for l in range(max(n,k), L+1):
                d = d + c_cos(n,k,l) * e2 ** l * varrho ** k * sp.cos(2 * n * psi)
    return d

def cos_phi_in_sin_pow(N: int, K: int) -> sp.core.Expr:
    """Symbolic sin-power series expansion of cos(phi)/cos(psi) - 1 up to given order.

    :param N: Maximum sin²-power.
    :param K: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for n in range(0, N+1):
        for k in range(1, K+1):
            for l in range(max(k,n), n+k):
                d = d + d_cos(n, k, l) * e2 ** l * varrho ** k * sin_psi ** (2 * n)
    return d

def sigma(J: int, delta: sp.core.Expr) -> sp.core.Expr:
    """ Symbolic "sigma" series up to given order.

    :param J: Maximum delta power.
    :param delta: Expression for delta.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for j in range(0, J, 2):
        j2 = j / 2
        d = d + sp.Rational(sp.Integer(-1) ** j2, sp.factorial(j)) * delta ** j2
    return d

def tau(J: int, omega: sp.core.Expr, delta: sp.core.Expr):
    """ Symbolic "tau" series up to given order.

    :param J: Maximum delta power.
    :param omega: Expression for "omega".
    :param delta: Expression for delta.
    :return: Symbolic representation of the series.
    """
    d = sp.Integer(0)
    for j in range(1, J, 2):
        j2 = (j - 1) / 2
        d = d + sp.Rational(sp.Integer(-1) ** j2, sp.factorial(j)) * delta ** j2
    return omega * d

def sin_phi_in_sin_pow2(N: int, K: int, J: int):
    """Symbolic sin-power series expansion of sin(phi)/sin(psi) up to given order.

    Using sigma and tau series.

    :param N: Maximum sin²-power.
    :param K: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    sin_psi2 = sin_psi ** 2
    o = phi_in_sin_pow(N, K)
    d = sin_psi2 * (sp.Integer(1) - sin_psi2) * o ** 2
    s = sigma(J, d)
    t = tau(J, o, d)
    return s + (sp.Integer(1) - sin_psi2) * t

def cos_phi_in_sin_pow2(N: int, K: int, J: int):
    """Symbolic sin-power series expansion of sin(phi)/sin(psi) up to given order.

    Using sigma and tau series.

    :param N: Maximum sin²-power.
    :param K: Maximum varrho power.
    :return: Symbolic representation of the series.
    """
    sin_psi2 = sin_psi ** 2
    o = phi_in_sin_pow(N, K)
    d = sin_psi2 * (sp.Integer(1) - sin_psi2) * o ** 2
    s = sigma(J, d)
    t = tau(J, o, d)
    return s - sin_psi2 * t

def phi_evo_sin_pow_dense(N, K):
    """Series for phi - pi/2 in sin powers for small rho with dense coefficients

    Dense in the meaning that there are no zero coefficients but instead the
    expressions for thw powers are more complicated.

    :param N: sin power limit.
    :param K: rho power limit.
    :return: symbolic series.
    """
    d = sp.S.Zero
    for n in range(0, N+1):
        for k in range(0, K+1):
            for l in range(0, n//2+k+1):
                d += cos_psi * sin_psi ** n * rho_ae2 ** (n + 1 + 2*k) * b_a ** ((n % 2) + 1 + 2*l) * d_phi_evo_dense(n, k, l)
    return d


def phi_evo_sin_pow(N, K):
    """Series for phi - pi/2 in sin-powers for small rho with simple sums

    The simple sums and exponents come at the cost of half of the computed
    coefficients being zero.

    :param N: sin power limit.
    :param K: rho powers limit.
    :return: Symbolic series.
    """
    d = sp.S.Zero
    for n in range(0, N+1):
        for k in range(n+1,K+1):
            for l in range(1,k+1):
                d += cos_psi * sin_psi ** n * rho_ae2 ** k * b_a ** l * c_phi_evo(n, k, l)
    return d


def phi_pow_evo_sin_pow(i, N, K):
    """Series for (phi - pi/2)^i in sin powers for small rho.

    :param i: integer power.
    :param N: sin power limit.
    :param K: rho power limit.
    :return: Symbolic series.
    """
    d = sp.S.Zero
    for n in range(0, N+1):
        for k in range(n+i, K+1):
            for r in range(i,k+1):
                d += sin_psi ** n * rho_ae2 ** k * b_a ** r * d_phi_pow_evo(n, k, r, i)
    return cos_psi ** i * d
