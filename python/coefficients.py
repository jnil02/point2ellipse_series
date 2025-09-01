
"""
Point-to-ellipse Fourier and sin-power series expansion coefficients.
"""

# External includes.
from math import floor, ceil
import sympy as sp

# Internal includes.
import series_substitutions
import polynomials
import cache
import symbols
from util import rf_half


@cache.ints_cache
def d_phi(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """(phi - psi) sin-power series expansion coefficients.

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
                    for t in range(l - n - k + m + r - p - q, min(l - k, r - q, l - n - k + m + r - q) + 1):
                        d = (d + sp.Rational(rf_half(k, r - q) * ((-1) ** (n - l - k) * 2 ** (r - 2 * q)),
                                             sp.factorial(q) * sp.factorial(r - 2 * q) * (m + 1 + r))
                             * sp.binomial(r - q, t) * sp.binomial(k - 1, m + r) * sp.binomial(m + 1, 2 * p + 1)
                             * sp.binomial(sp.Rational(k, 2) + r - q + l - k - t - 1, l - k - t)
                             * sp.binomial(p, n - (m + r - q + l - k - t - p)))
    return d

@cache.ints_cache
def d_phi2(n: int, k: int, l: int) -> sp.core.numbers.Rational:
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
                    for t in range(l - n - k + m + r - p - q, min(l - k, r - q) + 1):
                        d = (d + sp.Rational(rf_half(k, r - q) * ((-1) ** (n - l - k) * 2 ** (r - 2 * q)),
                                             sp.factorial(q) * sp.factorial(r - 2 * q) * (m + 1 + r))
                             * sp.binomial(r - q, t)
                             * sp.binomial(k - 1, m + r)
                             * sp.binomial(m + 1, 2 * p + 1)
                             * sp.binomial(sp.Rational(k, 2) + r - q + l - k - t - 1, l - k - t)
                             * sp.binomial(p + sp.Rational(1,2), n - (m + r - q + l - k - t - p)))
    return d

@cache.ints_cache
def c_phi(n: int, k: int, l: int) -> sp.core.numbers.Rational:
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
                            sum = sum + sp.binomial(2 * p + 1, i) * sp.binomial(1 + 2 * w, j) * sp.Integer((-1)**j)
                        for i in range(p-w+n, p + 1): # Empty if p-w+n > p.
                            j = w-p+i-n
                            sum = sum + sp.binomial(2 * p + 1, i) * sp.binomial(1 + 2 * w, j) * sp.Integer((-1)**j)
                        for i in range(p-w-n, p-n + 1):  # Empty if p-n < 0
                            j = w-p+i+n
                            sum = sum - sp.binomial(2 * p + 1, i) * sp.binomial(1 + 2 * w, j) * sp.Integer((-1)**j)
                        h = (h + sum * sp.Integer((-1) ** (l - k)) * sp.Rational(
                            rf_half(k, r - q),
                            sp.factorial(q) * sp.factorial(r - 2 * q) * (m + 1 + r) * 2 ** (2 * (m + l - k - t) + r + 1))
                            * sp.binomial(r - q, t)
                            * sp.binomial(k - 1, m + r) * sp.binomial(m + 1, 2 * p + 1)
                            * sp.binomial(sp.Rational(k, 2) + r - q + l - k - t - 1, l - k - t))
    return h

@cache.ints_cache
def d_phi_pow_polynomial(n: int, k: int, i: int) -> sp.core.Expr:
    # Polynomial for b_{n,i} in terms of {a_0,...,a_n}.
    tmp = series_substitutions.double_series_power_coeff(n, i)[k]
    # Polynomial for the varrho^k coefficient in b_{n,i} in terms of {a_{n,1},...a_{n,k+1}}
    tmp = series_substitutions.a_nk_sub(tmp, 1, d_phi)
    return tmp

@cache.ints_cache
def d_phi_pow(n: int, k: int, l: int, i: int) -> sp.core.numbers.Rational:
    """Compute (phi-psi)^i sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :param i: The power exponent.
    :return: Coefficient as a sympy rational number.
    """
    tmp = d_phi_pow_polynomial(n, k, i)
    return sp.expand(tmp).coeff(symbols.e2, l)  # Extract the l:th power of the series.

@cache.ints_cache
def d_sin_pow_polynomial(n: int, k: int, i: int) -> sp.core.Expr:
    # Polynomial for b_{n,i} in terms of {a_0,...,a_n}.
    tmp = series_substitutions.double_series_power_coeff(n, i)[k]
    # Polynomial for the delta^k coefficient in b_{n,i} in terms of {a_{n,1},...a_{n,k+1}}
    tmp = series_substitutions.a_nk_sub(tmp, 0, d_sin)
    return tmp

@cache.ints_cache
def d_sin_pow(n: int, k: int, l: int, i: int) -> sp.core.numbers.Rational:
    """Compute (sin(phi)/sin(psi)-1)^i sin-power series expansion coefficients.

    :param n: sin-power.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :param i: The power exponent.
    :return: Coefficient as a sympy rational number.
    """
    tmp = d_sin_pow_polynomial(n, k, i)
    return sp.expand(tmp).coeff(symbols.e2, l)  # Extract the l:th power of the series.

@cache.ints_cache
def d_N_nkl(n: int, k: int, l: int) -> sp.core.numbers.Rational:
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

@cache.ints_cache
def bp_nkl(n: int, k: int, l: int) -> sp.core.numbers.Rational:
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

@cache.ints_cache
def d_sin(n: int, k: int, l: int) -> sp.core.numbers.Rational:
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

@cache.ints_cache
def c_sin(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Compute sin(phi)/sin(psi)-1 Fourier series expansion coefficients.

    :param n: Fourier sin-multiple.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    return polynomials.sin_pow_to_cos_mul(n, k, l, 0, 0, d_sin)

@cache.ints_cache
def d_cos(n: int, k: int, l: int) -> sp.core.numbers.Rational:
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

@cache.ints_cache
def c_cos(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Compute cos(phi)/cos(psi)-1 Fourier series expansion coefficients.

    :param n: Fourier sin-multiple.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    return polynomials.sin_pow_to_cos_mul(n, k, l, 0, -1, d_cos)

@cache.ints_cache
def d_h(n: int, k: int, l: int) -> sp.core.numbers.Rational:
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

@cache.ints_cache
def c_h(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Compute (h+a-rho)/a Fourier series expansion coefficients.

    :param n: Fourier sin-multiple.
    :param k: varrho power of inner power series.
    :param l: e² power of innermost power series.
    :return: Coefficient as a sympy rational number.
    """
    return polynomials.sin_pow_to_cos_mul(n, k, l, 1, 0, d_h)
