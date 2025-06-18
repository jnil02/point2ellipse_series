
"""Point-to-ellipse Fourier and sin-power expansion coefficients.

For further details see
Nilsson, John-Olof. Point-to-ellipse Fourier series. DOI:
https://doi.org/10.48550/arXiv.2506.XXXXX.
"""

# External includes.
from math import floor, ceil
import sympy as sp
from functools import lru_cache

# Internal includes.
import polynomials
import util

@lru_cache(maxsize=None)
def d_phi(n, k, u):
    """Compute phi-psi sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    d = sp.Integer(0)
    for r in range(k - 1 + 1):
        for m in range(k - 1 - r + 1):
            for q in range(floor(r / 2) + 1):
                for p in range(floor(m / 2) + 1):
                    for t in range(min(u - k, r - q) + 1):
                        d = (d + sp.Rational(sp.Integer((-1) ** (2 * p + n - u - k)) * 2 ** (r - 2 * q) *
                            sp.functions.combinatorial.factorials.RisingFactorial(sp.Rational(k, 2), r - q),
                            sp.factorial(q) * sp.factorial(r - 2 * q) * (m + 1 + r))
                               * sp.binomial(r - q, t) * sp.binomial(k - 1, m + r) * sp.binomial(m + 1, 2 * p + 1)
                               * sp.binomial(sp.Rational(k, 2) + r - q + u - k - t - 1, u - k - t)
                               * sp.binomial(p, n - (m + r - q + u - k - t - p)))
    return d

@lru_cache(maxsize=None)
def d_phi2(n, k, u):
    """Compute phi-psi sin-power series expansion coefficients __with cos-sin factor integrated in the series__.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    d = sp.Integer(0)
    for r in range(k - 1 + 1):
        for m in range(k - 1 - r + 1):
            for q in range(floor(r / 2) + 1):
                for p in range(floor(m / 2) + 1):
                    for t in range(max(u - k, r - q) + 1):
                        d = (d + sp.Rational(sp.Integer((-1) ** (2*p+n-u - k)) * 2 ** (r - 2*q) *
                            sp.functions.combinatorial.factorials.RisingFactorial(sp.Rational(k, 2), r - q),
                            sp.factorial(q) * sp.factorial(r - 2 * q) * (m + 1 + r))
                               * sp.binomial(r - q, t) * sp.binomial(k - 1, m + r) * sp.binomial(m + 1, 2 * p + 1)
                               * sp.binomial(sp.Rational(k, 2) + r - q + u-k-t - 1, u-k-t)
                               * sp.binomial(p+sp.Rational(1,2), n - (m+r-q+u-k-t-p)))
    return d

@lru_cache(maxsize=None)
def c_phi(n, k, l):
    """Compute phi-psi Fourier series expansion coefficients.

    :param n: Fourier sin-multiple.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    h = sp.Integer(0)
    for r in range(k - 1 + 1):
        for m in range(k - 1 - r + 1):
            for p in range(floor(m / 2) + 1):
                for q in range(floor(r / 2) + 1):
                    for t in range(min(l - k, r - q) + 1):
                        w = m + r - q + l - k - t - p
                        sum = sp.Integer(0)
                        for i in range(max(0,p-n+1), min(p, p-n+1+w) + 1):  # Empty if w < 0.
                            j = w+p-i+1-n
                            sum = sum + sp.binomial(2*p+1, i) * sp.binomial(1 + 2*w, j) * sp.Integer((-1)**j)
                        for i in range(p-w+n, p + 1): # Empty if p-w+n > p.
                            j = w-p+i-n
                            sum = sum + sp.binomial(2*p+1, i) * sp.binomial(1 + 2*w, j) * sp.Integer((-1)**j)
                        for i in range(p-w-n, p-n + 1):  # Empty if p-n < 0
                            j = w-p+i+n
                            sum = sum - sp.binomial(2*p+1, i) * sp.binomial(1 + 2*w, j) * sp.Integer((-1)**j)
                        h = (h + sum * sp.Integer((-1) ** (l - k)) * sp.Rational(
                            sp.functions.combinatorial.factorials.RisingFactorial(sp.Rational(k, 2), r - q),
                            sp.factorial(q) * sp.factorial(r - 2 * q) * (m + 1 + r) * 2 ** (2 * (m + l - k - t) + r + 1))
                            * sp.binomial(r - q, t)
                            * sp.binomial(k - 1, m + r) * sp.binomial(m + 1, 2 * p + 1)
                            * sp.binomial(sp.Rational(k, 2) + r - q + l - k - t - 1, l - k - t))
    return h

@lru_cache(maxsize=None)
def d_phi_pow(n, k, l, i):
    """Compute (phi-psi)^i sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :param i: The power exponent.
    :return: Coefficient as a sympy rational number.
    """
    # NOTE, identical to d_sin_pow except a_nk_sub arguments.
    # Polynomial for b_{n,i} in terms of {a_0,...,a_n}.
    b_ni = polynomials.power_of_power_series_coefficient_polynomial(n, i, 'a_')
    # Polynomial for the varrho^k coefficient in b_{n,i} in terms of {a_{n,1},...a_{n,k+1}}
    tmp = util.poly_bell_substitution(b_ni, 'a')[k]
    tmp = util.a_nk_sub(tmp, 1, d_phi)
    return sp.expand(tmp).coeff(util.e2, l)  # Extract the l:th power of the series.

@lru_cache(maxsize=None)
def d_sin_pow(n, k, l, i):
    """Compute (sin(phi)/sin(psi)-1)^i sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :param i: The power exponent.
    :return: Coefficient as a sympy rational number.
    """
    # NOTE, identical to d_phi_pow except a_nk_sub arguments.
    # Polynomial for b_{n,i} in terms of {a_0,...,a_n}.
    b_ni = polynomials.power_of_power_series_coefficient_polynomial(n, i, 'a_')
    # Polynomial for the delta^k coefficient in b_{n,i} in terms of {a_{n,1},...a_{n,k+1}}
    tmp = util.poly_bell_substitution(b_ni, 'a')[k]
    tmp = util.a_nk_sub(tmp, 0, d_sin)
    return sp.expand(tmp).coeff(util.e2, l)  # Extract the l:th power of the series.

@lru_cache(maxsize=None)
def d_N_nkl(n,k,l):
    """Compute inverse radius of curvature sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    d = sp.Integer(0)
    for i in range(1,min(n,l)+1):
        for j in range(min(2*i,k)+1):
            d = (d + sp.binomial(sp.Rational(1,2), i)
                 * sp.binomial(2*i,j) * (-1) ** i * d_sin_pow(n - i, k, l - i, j))
    return d

@lru_cache(maxsize=None)
def bp_nkl(n,k,l):
    """Compute cos(phi-psi) sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    c = sp.Integer(0)
    for i in range(1, min(n,floor(k/2))+1):
        for j in range(min(i,n-i)+1):
            c = c + (-1) ** (i+j) * sp.binomial(i,j) / sp.factorial(2*i) * d_phi_pow(n - i - j, k, l, 2 * i)
    return c

@lru_cache(maxsize=None)
def d_sin(n, k, l):
    """Compute sin(phi)/sin(psi)-1 sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    d = sp.Integer(0)
    for i in range(1,min(k,2*n+1)+1):
        for j in range(max(0,ceil(i/2)-l+n), min(ceil(i/2),n-floor(i/2))+1):
            d = d + sp.Rational(sp.binomial(ceil(i/2), j) * (-1) ** (floor(i/2)+j)
                                * d_phi_pow(n - floor(i / 2) - j, k, l, i), sp.factorial(i))
    return d

@lru_cache(maxsize=None)
def c_sin(n,k,l):
    """Compute sin(phi)/sin(psi)-1 Fourier series expansion coefficients.

    :param n: Fourier sin-multiple.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    return util.sin_pow_to_cos_mul(n, k, l, d_sin, 0, 0)

@lru_cache(maxsize=None)
def d_cos(n, k, l):
    """Compute cos(phi)/cos(psi)-1 sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    d = sp.Integer(0)
    for i in range(1,min(k,2*n+1)+1):
        for j in range(max(0,floor(i/2)-l+n), min(floor(i/2),n-ceil(i/2))+1):
            d = d + sp.Rational(sp.binomial(floor(i/2), j) * (-1) ** (ceil(i/2)+j)
                                * d_phi_pow(n - ceil(i / 2) - j, k, l, i), sp.factorial(i))
    return d

@lru_cache(maxsize=None)
def c_cos(n,k,l):
    """Compute cos(phi)/cos(psi)-1 Fourier series expansion coefficients.

    :param n: Fourier sin-multiple.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    return util.sin_pow_to_cos_mul(n, k, l, d_cos, 0, -1)

@lru_cache(maxsize=None)
def d_h(n, k, l):
    """Compute (h+a-rho)/a sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    if k==0:
        return -d_N_nkl(n,k,l)
    else:
        return bp_nkl(n,k+1,l) - d_N_nkl(n,k,l)

@lru_cache(maxsize=None)
def c_h(n,k,l):
    """Compute (h+a-rho)/a Fourier series expansion coefficients.

    :param n: Fourier sin-multiple.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    return util.sin_pow_to_cos_mul(n, k, l, d_h, 1, 0)
