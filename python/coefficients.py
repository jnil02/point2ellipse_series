"""
Point-to-ellipse far-field series expansion coefficients.
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
    """Coefficient of the sin-power series for the normal direction φ − ψ.

    Far-field (exterior) base coefficient, from which the sin(φ), cos(φ) and
    height series far from the center are built.

    :param n: sin(ψ) power.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 0 and k >= 1 and l >= max(n + 1, k) and l <= n + k, f"d_phi indices out of range. n: {n} k: {k} l: {l}"
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
    """Coefficient of the sin-power series for φ − ψ, cos·sin factor folded in.

    Variant of d_phi with the cos·sin prefactor absorbed into the series.

    :param n: sin(ψ) power.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 0 and k >= 1 and l >= k and l <= n + k, f"d_phi2 indices out of range. n: {n} k: {k} l: {l}"
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
                             * sp.binomial(p + sp.Rational(1, 2), n - (m + r - q + l - k - t - p)))
    return d


@cache.ints_cache
def c_phi(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Coefficient of the Fourier series for the normal direction φ − ψ.

    :param n: Fourier multiple index.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 1 and k >= 1 and l >= max(n, k), f"c_phi indices out of range. n: {n} k: {k} l: {l}"
    h = sp.Integer(0)
    for r in range(k - 1 + 1):
        for m in range(k - 1 - r + 1):
            for p in range(floor(m / 2) + 1):
                for q in range(floor(r / 2) + 1):
                    for t in range(min(l - k, r - q) + 1):
                        w = m + r - q + l - k - t - p
                        sum = sp.Integer(0)
                        for i in range(max(0, p - n + 1), min(p, p - n + 1 + w) + 1):  # Empty if w < 0.
                            j = w + p - i + 1 - n
                            sum = sum + sp.binomial(2 * p + 1, i) * sp.binomial(1 + 2 * w, j) * sp.Integer((-1) ** j)
                        for i in range(p - w + n, p + 1):  # Empty if p-w+n > p.
                            j = w - p + i - n
                            sum = sum + sp.binomial(2 * p + 1, i) * sp.binomial(1 + 2 * w, j) * sp.Integer((-1) ** j)
                        for i in range(p - w - n, p - n + 1):  # Empty if p-n < 0
                            j = w - p + i + n
                            sum = sum - sp.binomial(2 * p + 1, i) * sp.binomial(1 + 2 * w, j) * sp.Integer((-1) ** j)
                        h = (h + sum * sp.Integer((-1) ** (l - k)) * sp.Rational(
                            rf_half(k, r - q),
                            sp.factorial(q) * sp.factorial(r - 2 * q) * (m + 1 + r) * 2 ** (
                                        2 * (m + l - k - t) + r + 1))
                             * sp.binomial(r - q, t)
                             * sp.binomial(k - 1, m + r) * sp.binomial(m + 1, 2 * p + 1)
                             * sp.binomial(sp.Rational(k, 2) + r - q + l - k - t - 1, l - k - t))
    return h


@cache.ints_cache
def d_phi_pow_polynomial(n: int, k: int, i: int) -> sp.core.Expr:
    # Polynomial for b_{n,i} in terms of {a_0,...,a_n}.
    tmp = series_substitutions.double_series_power_coeff(n, i)[k]
    # Polynomial for the varrho^k coefficient in b_{n,i} in terms of {a_{n,1},...a_{n,k+1}}
    tmp = series_substitutions.a_nk_sub(tmp, lambda n, k: series_substitutions.a_nk_ser(n, k, 1, d_phi, symbols.e2))
    return tmp


@cache.ints_cache
def d_phi_pow(n: int, k: int, l: int, i: int) -> sp.core.numbers.Rational:
    """Coefficient of the sin-power series for (φ − ψ)^i.

    Powers of the φ − ψ series, used to build the sin(φ), cos(φ) and
    sin(φ)/sin(ψ) series via their Taylor expansions.

    :param n: sin(ψ) power.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :param i: power exponent.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 0 and k >= i and l >= max(n + i, k) and l <= n + k and i >= 1, \
        f"d_phi_pow indices out of range. n: {n} k: {k} l: {l} i: {i}"
    tmp = d_phi_pow_polynomial(n, k, i)
    return sp.expand(tmp).coeff(symbols.e2, l)  # Extract the l:th power of the series.


@cache.ints_cache
def d_sin_pow_polynomial(n: int, k: int, i: int) -> sp.core.Expr:
    # Polynomial for b_{n,i} in terms of {a_0,...,a_n}.
    tmp = series_substitutions.double_series_power_coeff(n, i)[k]
    # Polynomial for the delta^k coefficient in b_{n,i} in terms of {a_{n,1},...a_{n,k+1}}
    tmp = series_substitutions.a_nk_sub(tmp, lambda n, k: series_substitutions.a_nk_ser(n, k, 0, d_sin, symbols.e2))
    return tmp


@cache.ints_cache
def d_sin_pow(n: int, k: int, l: int, i: int) -> sp.core.numbers.Rational:
    """Coefficient of the sin-power series for (sin(φ)/sin(ψ) − 1)^i.

    :param n: sin(ψ) power.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :param i: power exponent.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 0 and k >= 0 and l >= 0 and l <= n + k and i >= 0, \
        f"d_sin_pow indices out of range. n: {n} k: {k} l: {l} i: {i}"
    tmp = d_sin_pow_polynomial(n, k, i)
    return sp.expand(tmp).coeff(symbols.e2, l)  # Extract the l:th power of the series.


@cache.ints_cache
def d_N_nkl(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Coefficient of the sin-power series for the inverse radius of curvature.

    :param n: sin(ψ) power.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 1 and k >= 0 and l >= max(n, k + 1) and l <= n + k, \
        f"d_N_nkl indices out of range. n: {n} k: {k} l: {l}"
    d = sp.Integer(0)
    for i in range(1, n + 1):
        for j in range(min(2 * i, k) + 1):
            d = (d + sp.binomial(sp.Rational(1, 2), i)
                 * sp.binomial(2 * i, j) * (-1) ** i * d_sin_pow(n - i, k, l - i, j))
    return d


@cache.ints_cache
def bp_nkl(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Coefficient of the sin-power series for cos(φ − ψ).

    :param n: sin(ψ) power.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 1 and k >= 1 and l >= max(n, k) and l <= n + k - 1, \
        f"bp_nkl indices out of range. n: {n} k: {k} l: {l}"
    c = sp.Integer(0)
    for i in range(1, min(n, floor(k / 2)) + 1):
        for j in range(max(0, n - l + i), min(i, n - i, n + k - l - i) + 1):
            c = c + (-1) ** (i + j) * sp.binomial(i, j) / sp.factorial(2 * i) * d_phi_pow(n - i - j, k, l, 2 * i)
    return c


@cache.ints_cache
def d_sin(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Coefficient of the sin-power series for sin(φ)/sin(ψ) − 1.

    :param n: sin(ψ) power.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 0 and k >= 1 and l >= max(n, k) and l <= n + k, \
        f"d_sin indices out of range. n: {n} k: {k} l: {l}"
    d = sp.Integer(0)
    for i in range(1, min(k, 2 * n + 1) + 1):
        for j in range(max(0, ceil(i / 2) - l + n), min(ceil(i / 2), n - floor(i / 2), n + k - l - floor(i / 2)) + 1):
            d = d + sp.Rational(sp.binomial(ceil(i / 2), j) * (-1) ** (floor(i / 2) + j)
                                * d_phi_pow(n - floor(i / 2) - j, k, l, i), sp.factorial(i))
    return d


@cache.ints_cache
def c_sin(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Coefficient of the Fourier series for sin(φ)/sin(ψ) − 1.

    :param n: Fourier multiple index.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 0 and k >= 1 and l >= max(n, k), \
        f"c_sin indices out of range. n: {n} k: {k} l: {l}"
    return polynomials.sin_pow_to_cos_mul(n, k, l, 0, 0, d_sin)


@cache.ints_cache
def d_cos(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Coefficient of the sin-power series for cos(φ)/cos(ψ) − 1.

    :param n: sin(ψ) power.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 0 and k >= 1 and l >= max(n, k) and l <= n + k - 1, \
        f"d_cos indices out of range. n: {n} k: {k} l: {l}"
    d = sp.Integer(0)
    for i in range(1, min(k, 2 * n + 1) + 1):
        for j in range(max(0, floor(i / 2) - l + n), min(floor(i / 2), n - ceil(i / 2), n + k - l - ceil(i / 2)) + 1):
            d = d + sp.Rational(sp.binomial(floor(i / 2), j) * (-1) ** (ceil(i / 2) + j)
                                * d_phi_pow(n - ceil(i / 2) - j, k, l, i), sp.factorial(i))
    return d


@cache.ints_cache
def c_cos(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Coefficient of the Fourier series for cos(φ)/cos(ψ) − 1.

    :param n: Fourier multiple index.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 0 and k >= 1 and l >= max(n, k), \
        f"c_cos indices out of range. n: {n} k: {k} l: {l}"
    return polynomials.sin_pow_to_cos_mul(n, k, l, 0, -1, d_cos)


@cache.ints_cache
def d_h(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Coefficient of the sin-power series for the height term (h + a − ρ)/a.

    :param n: sin(ψ) power.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 1 and k >= 0 and l >= max(n, k + 1) and l <= n + k, \
        f"d_h indices out of range. n: {n} k: {k} l: {l}"
    if k == 0:
        return -d_N_nkl(n, k, l)
    else:
        return bp_nkl(n, k + 1, l) - d_N_nkl(n, k, l)


@cache.ints_cache
def c_h(n: int, k: int, l: int) -> sp.core.numbers.Rational:
    """Coefficient of the Fourier series for the height term (h + a − ρ)/a.

    :param n: Fourier multiple index.
    :param k: ϱ (= a/ρ) power.
    :param l: e² power.
    :return: Coefficient as a sympy rational.
    """
    assert n >= 0 and k >= 0 and l >= max(n, k + 1), \
        f"c_h indices out of range. n: {n} k: {k} l: {l}"
    return polynomials.sin_pow_to_cos_mul(n, k, l, 1, 0, d_h)
