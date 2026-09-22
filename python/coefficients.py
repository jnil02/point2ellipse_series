"""
Point-to-ellipse Fourier and sin-power series expansion coefficients.
"""

# External includes.
from math import floor, ceil
import sympy as sp
from sympy import rf
from sympy.functions.combinatorial.numbers import stirling

# Internal includes.
import series_substitutions
import polynomials
import cache
import symbols
from util import rf_half, E2


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
def d_phi_evo(k: int, l: int, n: int) -> sp.core.numbers.Rational:
    """Dense φ−π/2 coefficient (single series in σ).

    Coefficient of the series for the ellipse-normal direction φ (the geodetic
    latitude) at a point inside the evolute. Dense form of c_phi_evo (vanishing
    terms removed).

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    s = k % 2
    m_k = (k - 1) // 2
    d = sp.S.Zero
    for j in range(n + 1):
        for q in range(2 * j, 1 - s + 2 * n + 1):
            d += 2 ** (q - 2 * j) * sp.binomial(q - j, j) * sp.binomial(k - 2 - q, 1 - s + 2 * n - q) * sp.binomial(
                m_k - n, l - n + j)
    return sp.Rational(d * (-1) ** (l + n + k), k) * sp.binomial(sp.Rational(k, 2), m_k - n)


@cache.ints_cache
def c_phi_evo(k, l, n):
    """Sparse φ−π/2 coefficient (every other term vanishes by parity).

    Coefficient of the series for the ellipse-normal direction φ (the geodetic
    latitude) at a point inside the evolute; the base coefficient from which the
    sin(φ), cos(φ) and h series are built.

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= l + 1 and n >= 1 and n <= k, \
        f"c_phi_evo indices out of range. k: {k} l: {l} n: {n}"
    c = sp.S.Zero
    if k < 1:
        return c
    if (k - l - 1) % 2 != 0:
        return c
    if (n - k) % 2 != 0:
        return c
    for j in range(floor((n - 1) / 2) + 1):
        b = sp.binomial((k - n) // 2, (l - n + 1) // 2 + j)
        s = sp.S.Zero
        for q in range(2 * j, n - 1 + 1):
            # Note, the first binomial can become (-1,0) which should be 1.
            s += 2 ** (q - 2 * j) * sp.binomial(k - 2 - q, n - 1 - q) * sp.binomial(q - j, j)
        c += sp.binomial(sp.Rational(k, 2), (k - n) // 2) * s * b
    return c * (-1) ** ((l + n + 1) // 2) / k


@cache.ints_cache
def d_phi_pow_evo_polynomial(k: int, l: int, i: int) -> sp.core.Expr:
    # Polynomial for A_{l,i} in terms of {a_0,...,a_l}. The inside-evolute
    # generator a_l starts at z^{l+1}, so use the tightened variant based on
    # the nonzero condition k>=i*(l+1)).
    tmp = series_substitutions.double_series_power_coeff_evo(l, i)[k]
    # Polynomial for the rho^k coefficients in A_{l,i} in terms of {a_{l,1},...a_{l,k+1}}
    tmp = series_substitutions.a_nk_sub(tmp, lambda l, k: series_substitutions.a_nk_C(l, k, c_phi_evo, symbols.e2))
    return tmp


@cache.ints_cache
def c_phi_pow_evo(k: int, l: int, n: int, i: int) -> sp.core.numbers.Rational:
    """Sparse (φ−π/2)^i coefficient.

    Coefficient of the i-th power of the φ−π/2 series, used to build the
    sin(φ), cos(φ) and sin(ψ)/sin(φ) series via their Taylor expansions.

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :param i: power exponent.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= l + i and n >= i and n <= k, \
        f"c_phi_pow_evo indices out of range. k: {k} l: {l} n: {n} i: {i}"
    if (i + l - k) % 2 != 0 or (n - k) % 2 != 0:  # Parity constraint from the underlying coefficients.
        return sp.S.Zero
    tmp = d_phi_pow_evo_polynomial(k, l, i)
    return sp.expand(tmp).coeff(symbols.e2, n)  # Extract the n:th e2-power of the series.


@cache.ints_cache
def c_sin_phi_evo(k: int, l: int, n: int) -> sp.core.numbers.Rational:
    """Sparse sin(φ) coefficient.

    Coefficient of the series for sin(φ), the vertical component of the unit
    ellipse normal (the n-vector) at a point inside the evolute.

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= l and n >= 0 and n <= k, \
        f"c_sin_phi_evo indices out of range. k: {k} l: {l} n: {n}"
    if (l - k) % 2 != 0 or (l - n) % 2 != 0:
        return sp.S.Zero
    d = sp.S.Zero
    for i in range(n // 2 + 1):
        for j in range(max(0, ceil((l + 2 * i - k) / 2.)), min(i, l // 2) + 1):
            d += (sp.Rational((-1) ** (i + j), sp.factorial(2 * i)) * sp.binomial(i, j)
                  * c_phi_pow_evo(k, l - 2 * j, n, 2 * i))
    return d


@cache.ints_cache
def d_sin_phi_evo(k: int, l: int, n: int) -> sp.core.numbers.Rational:
    """Dense sin(φ) coefficient.

    Dense-form coefficient of the sin(φ) series (see c_sin_phi_evo).

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= 2 and n >= 1 and l <= k // 2 and n <= k // 2, \
        f"d_sin_phi_evo indices out of range. k: {k} l: {l} n: {n}"
    return c_sin_phi_evo(k, 2 * l + (k % 2), 2 * n + (k % 2))


@cache.ints_cache
def c_cos_phi_evo(k: int, l: int, n: int) -> sp.core.numbers.Rational:
    """Sparse cos(φ)/|cos(ψ)| coefficient.

    Coefficient of the series for cos(φ), the horizontal normal component, with
    the |cos(ψ)| factor extracted so the expansion stays regular at ψ = ±π/2.

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= l and n >= 1 and n <= k, \
        f"c_cos_phi_evo indices out of range. k: {k} l: {l} n: {n}"
    if (l + 1 - k) % 2 != 0 or (n - k) % 2 != 0:
        return sp.S.Zero
    d = sp.S.Zero
    for i in range((n - 1) // 2 + 1):
        for j in range(max(0, ceil((l + 2 * i + 1 - k) / 2.)), min(i, l // 2) + 1):
            d += (sp.Rational((-1) ** (i + 1 + j), sp.factorial(2 * i + 1)) * sp.binomial(i, j)
                  * c_phi_pow_evo(k, l - 2 * j, n, 2 * i + 1))
    return d


@cache.ints_cache
def d_cos_phi_evo(k: int, l: int, n: int) -> sp.core.numbers.Rational:
    """Dense cos(φ)/|cos(ψ)| coefficient.

    Dense-form coefficient of the cos(φ)/|cos(ψ)| series (see c_cos_phi_evo).

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert k >= 1 and l >= 0 and n >= (k + 1) % 2 and l <= (k - 1) // 2 and n <= ceil((k - 1) / 2.), \
        f"d_cos_phi_evo indices out of range. k: {k} l: {l} n: {n}"
    p = (k - 1) % 2
    return c_cos_phi_evo(k, p + 2 * l, 2 * n + 1 - p)


@cache.ints_cache
def c_sin_phi_inv_evo(k: int, l: int, n: int) -> sp.core.numbers.Rational:
    """Sparse sin(ψ)/sin(φ) coefficient.

    Coefficient of the intermediate sin(ψ)/sin(φ) series used in deriving the
    height h.

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= l and n >= 0 and n <= k, \
        f"c_sin_phi_inv_evo indices out of range. k: {k} l: {l} n: {n}"
    if (l - k) % 2 != 0 or (l - n) % 2 != 0:
        return sp.S.Zero
    d = sp.S.Zero
    for i in range(n // 2 + 1):
        for j in range(max(0, ceil((l + 2 * i - k) / 2.)), l // 2 + 1):
            d += (sp.Rational(E2(i) * (-1) ** (j), sp.factorial(2 * i)) * sp.binomial(i, j)
                  * c_phi_pow_evo(k, l - 2 * j, n, 2 * i))
    return d


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


@cache.ints_cache
def R(k: int, l: int, n: int, i: int) -> sp.core.numbers.Rational:
    """Intermediate alternating sum over (φ−π/2)^{2i} coefficients feeding c_N_evo.

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :param i: power exponent.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= l and n >= 0 and n <= k and i >= 0 and i <= n // 2, \
        f"R indices out of range. k: {k} l: {l} n: {n} i: {i}"
    s = sp.S.Zero
    for j in range(max(0, ceil((l + 2 * i - k) / 2.)), l // 2 + 1):
        s += (-1) ** (j) * sp.binomial(i, j) * c_phi_pow_evo(k, l - 2 * j, n, 2 * i)
    return s


@cache.ints_cache
def c_N_evo(k: int, l: int, n: int) -> sp.core.Rational:
    """Sparse ε²·N/a coefficient (N the normal radius of curvature).

    Coefficient of the series for the radius-of-curvature term that enters the
    height h.

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= l and n >= 0 and n <= k + 1, \
        f"c_N_evo indices out of range. k: {k} l: {l} n: {n}"
    if (l - k) % 2 != 0 or (l - n - 1) % 2 != 0:
        return sp.S.Zero
    d = sp.S.Zero
    for p in range(k + 1):
        for i in range(p // 2 + 1):
            t = p + 1 - n
            if (t % 2 == 0 and t >= 0 and t <= 2 * i):
                d += C_mt(i, t // 2) * R(k, l, p, i)
    return d


@cache.ints_cache
def cp_evo(k: int, l: int, n: int) -> sp.core.Rational:
    """Sparse sin(ψ)/sin(φ) contribution to h, with ρ·sin(ψ)/a pulled out (feeds c_h_evo).

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 1 and k >= l and n >= 1 and n <= k + 1, \
        f"cp_evo indices out of range. k: {k} l: {l} n: {n}"
    if (l - k) % 2 != 0 or (k + 1 - n) % 2 != 0:
        return sp.S.Zero
    if n <= 2 and n <= k - 1:
        # k = 0, l=0, n = 2
        return c_sin_phi_inv_evo(k - 1, l - 1, n)
    if 3 <= n and n <= k - 1:
        return c_sin_phi_inv_evo(k - 1, l - 1, n) - c_sin_phi_inv_evo(k - 1, l - 1, n - 2)
    if 3 <= n and n <= k + 1:
        return -c_sin_phi_inv_evo(k - 1, l - 1, n - 2)
    # This will happen for e.g. 1,1,2.
    return sp.S.Zero


@cache.ints_cache
def c_h_evo(k: int, l: int, n: int) -> sp.core.Rational:
    """Sparse h/a coefficient with ρ·sin(ψ)/a pulled out (every other term vanishes).

    Coefficient of the series for the signed distance (height) h from the
    ellipse to the point; combines the curvature and sin(ψ)/sin(φ) parts
    (c_h_evo = cp_evo − c_N_evo).

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= l and n >= 1 and n <= k + 1, \
        f"c_h_evo indices out of range. k: {k} l: {l} n: {n}"
    if (l - k) % 2 != 0 or (l - n - 1) % 2 != 0:
        return sp.S.Zero
    if l == 0:
        return -c_N_evo(k, l, n)
    return cp_evo(k, l, n) - c_N_evo(k, l, n)


@cache.ints_cache
def d_h_evo(k: int, l: int, n: int) -> sp.core.numbers.Rational:
    """Dense h/a coefficient with ρ·sin(ψ)/a pulled out.

    Dense-form coefficient of the height (h) series (see c_h_evo).

    :param k: σ (= ρ/(a·e²)) power.
    :param l: sin(ψ) power.
    :param n: ε (= b/a) power.
    :return: Coefficient as a sympy rational.
    """
    assert l >= 0 and k >= 0 and n >= 0 and n <= k // 2 and l <= k // 2, \
        f"d_h_evo indices out of range. k: {k} l: {l} n: {n}"
    sn = k % 2
    return c_h_evo(k, 2 * l + sn, 2 * n + 1 + sn)


@cache.ints_cache
def a_mr(m: int, r: int) -> sp.core.numbers.Rational:
    """Helper coefficient for the ε²·N/a finite-sum reduction.

    Turns the infinite e^{2m} series in N into a finite polynomial in ε
    (together with B_rt and C_mt).

    :param m: outer index.
    :param r: inner index (0 <= r <= m).
    :return: Coefficient as a sympy rational.
    """
    assert r <= m and r >= 0, f"a_mr index out of range. r: {r} m: {m}"
    if m == 0 and r == 0:
        return sp.S.One
    a = sp.S.Zero
    for k in range(r, m + 1):
        b = sp.S.Zero
        for t in range(1, k + 1):
            b += (-1) ** (k - t) * sp.binomial(2 * k, k - t) * t ** (2 * m)
        a += sp.Rational(1, sp.factorial(k) * 2 ** (k - r)) * b * stirling(k, r, kind=1, signed=True)
    return sp.Rational(2 * (-1) ** m, sp.factorial(2 * m)) * a


@cache.ints_cache
def C_mt(m: int, t: int) -> sp.core.numbers.Rational:
    """Helper coefficient (combination of a_mr and B_rt) for the ε²·N/a reduction.

    :param m: outer index.
    :param t: inner index (0 <= t <= m).
    :return: Coefficient as a sympy rational.
    """
    assert m >= 0 and t >= 0 and t <= m, \
        f"C_mt indices out of range. m: {m} t: {t}"
    C = sp.S.Zero
    for r in range(t, m + 1):
        C += a_mr(m, r) * B_rt(r, t)
    return C


@cache.ints_cache
def B_rt(r: int, t: int) -> sp.core.numbers.Rational:
    """Helper coefficient (Stirling/rising-factorial sum) for the ε²·N/a reduction.

    :param r: outer index.
    :param t: inner index (0 <= t <= r).
    :return: Coefficient as a sympy rational.
    """
    assert r >= t and t >= 0, f"B_rt indices out of range. r: {r} t: {t}"
    B = sp.S.Zero
    for k in range(t, r + 1):
        B += stirling(r, k, kind=2) * rf(sp.Rational(1, 2), k) * (-1) ** (k - t) * sp.binomial(k, t)
    return B
