
"""Various substitutions of series into polynomials.
"""

# External includes.
from collections.abc import Callable
import sympy as sp

# Internal includes.
import polynomials
import cache
import series


def poly_bell_substitution(p: sp.core.Expr,
                           start: Callable[[int], int] = lambda l: 1) -> series.SeriesBase:
    r"""Expand a polynomial p(a_0,...,a_n) by substituting Bell polynomials for a_n^i.

    p is a multidimensional polynomial in variables a_n. (The assumption is
    that the symbol name ends with "_<number>" which can be parsed and
    transferred to the substitution variable.)
    Assume

    a_n = \sum_{k=start(n)}^\infty b_{n,k} z^k

    where b_{n,k} are new variables, then

    a_n^i = \sum_{k=i*start(n)} \hat{B}_{k,i}(b_{n,1},...,b_{n,k-i+1})z^k

    where \hat{B}_{k,i}(...) is a partial ordinary Bell polynomial.
    See https://en.wikipedia.org/wiki/Bell_polynomials

    This function makes the Bell polynomial substitution.

    :param p: The polynomial to make the substitution in. Polynomial must be expanded.
    :param start: Lowest power z^{start(n)} at which generator a_n begins. The
        default assumes every generator starts at z^1, giving the standard partial
        ordinary Bell nonzero condition k>=i. For generators that start at
        z^{n+1} (e.g. the inside-evolute a_l), pass start=lambda l: l+1, which
        tightens the nonzero condition to k>=i*(n+1). Since b_{n,k}=0 for
        k<start(n), the extra terms are identically zero, so this only prunes
        provably-zero terms and does not change the result.
    :return: A coefficient sequence representing the substituted polynomial as a series.
    """
    seqTot = series.SeriesEmpty()  # Sequence for the whole polynomial.
    # Split the polynomial in terms, handle each term and add the results.
    for polyTerm in sp.Add.make_args(sp.expand(p)):
        seqTerm = series.SeriesFactor(1)  # Sequence for term.
        # Split the term in factors, handle each factor and multiply the results.
        for termFactor in sp.Mul.make_args(polyTerm):
            # Handle each factor type.
            if termFactor.is_Number:
                seqTerm = seqTerm * termFactor
            elif termFactor.is_Pow or termFactor.is_Symbol:
                # Handling of x^y and x^1.
                base_exp = termFactor.as_base_exp() if termFactor.is_Pow else (termFactor, 1)
                # Retrieve the indices of the coefficient, e.g. "_0" for "a_0".
                ix = base_exp[0].name[base_exp[0].name.find('_'):]
                l_gen = int(ix[1:])  # Generator index l from the symbol name a_l.
                # The i:th power of a_l starts at z^{i*start(l)}; below that the
                # partial ordinary Bell polynomial vanishes, so the guard uses it.
                def bellLambdaGen(j, x, thr):
                    return lambda n: polynomials.partial_ordinary_bell_polynomial(n, j, x) if n >= thr else 0
                seqTerm = seqTerm * series.Series(
                    bellLambdaGen(int(base_exp[1]), 'a' + ix, int(base_exp[1]) * start(l_gen)))
            else:
                raise Exception("Unhandled factor in sympy expression:" + str(termFactor))
        seqTot = seqTot + seqTerm
    return seqTot

@cache.ints_cache
def double_series_power_coeff(n: int, i: int) -> series.SeriesBase:
    """Coefficient of a power of a double power series (outer index from 0, inner from 1).

    Far-field counterpart of double_series_power_coeff_evo.

    :param n: first index of the resulting series coefficients.
    :param i: power of the double power series.
    :return: Series representing the coefficient.
    """
    # Polynomial for b_{n,i} in terms of {a_0,...,a_n}.
    b_ni = polynomials.ordinary_potential_polynomial(n, i, "a")
    # Polynomial for the varrho^k coefficient in b_{n,i} in terms of {a_{n,1},...a_{n,k+1}}
    return poly_bell_substitution(b_ni)

@cache.ints_cache
def double_series_power_coeff_evo(l: int, i: int) -> series.SeriesBase:
    """Coefficient of the power of a double power series for the inside-evolute chain.

    Same as double_series_power_coeff, but the inner generator a_l starts at
    z^{l+1} rather than z^1. This tightens the Bell substitution's nonzero
    condition to k>=i*(l+1) (manuscript: the kappa-sum defining c^phi_{k,l,n,i}
    starts at i_l(l+1)). The result is unchanged from the untightened version;
    the pruned terms are identically zero after the a_{l,k} substitution.

    :param l: sin-power (first index of the resulting series coefficients).
    :param i: The power of the double power series.
    :return: Series representing the coefficient.
    """
    b_li = polynomials.ordinary_potential_polynomial(l, i, "a")
    return poly_bell_substitution(b_li, start=lambda l: l + 1)

def a_nk_ser(n: int, k: int, n_offset: int, d_nkl: Callable[[int, int, int], sp.core.Expr], e2: sp.core.Symbol) -> sp.core.Expr:
    """Far-field generator a_{n,k}: the finite e²-series sum_l d_nkl(n,k,l)*e2^l.

    Far-field counterpart of a_nk_C; l runs over max(k, n+n_offset) .. n+k.

    :param n: sin(ψ) power (generator index).
    :param k: ϱ power.
    :param n_offset: offset setting the lowest e² power l.
    :param d_nkl: coefficient callback d_nkl(n, k, l).
    :param e2: series base variable (e²).
    :return: Finite series expression.
    """
    a_nk = sp.S.Zero
    for l in range(max(k, n + n_offset), n + k + 1):
        a_nk = a_nk + d_nkl(n, k, l) * e2 ** l
    return a_nk

def a_nk_C(l: int, k: int, c: Callable[[int, int, int], sp.core.Expr], e2: sp.core.Symbol):
    """Specific finite a_{l,k} series from 1 to k if k>=l+1.

    The generator a_{l,k} carries the parity constraints of the underlying c
    coefficients: it is zero unless k >= l+1 and k == l+1 (mod 2), and within a
    nonzero generator only the e2-powers n == k (mod 2) contribute. Off-parity
    generators/terms are identically zero, so pruning them here does not change
    the result, only avoids evaluating provably-zero coefficients.

    :param l: sin-power index (the generator index).
    :param k: sigma-power index.
    :param c: coefficient callback c(k, l, n).
    :param e2: series base variable.
    :return: expression for finite series.
    """
    a_lk = sp.S.Zero
    if k < l + 1 or (k - l - 1) % 2 != 0:
        return a_lk
    for n in range(2 - k % 2, k + 1, 2):  # e2-power n == k (mod 2); others vanish
        a_lk += c(k, l, n) * e2 ** n
    return a_lk

def a_nk_sub(p: sp.core.Expr, a_nk: Callable[[int, int], sp.core.Expr]) -> sp.core.Expr:
    """Substitute a_{n,k} with the result of the callback.

    The indices n and k are carried over from the symbol names to the callback arguments.

    :param a_nk: Function for retrieving the a_nk coefficient substitutions.
    :param p: The polynomial to make the substitution in.
    :return: The polynomial with the substitution.
    """
    pTot = sp.Integer(0)
    # Split the polynomial in terms, handle each term and add the results.
    for polyTerm in sp.Add.make_args(sp.expand(p)):
        pTerm = 1
        # Split the term in factors, handle each factor and multiply the results.
        for termFactor in sp.Mul.make_args(polyTerm):
            # Handle each factor type.
            if termFactor.is_Number:
                pTerm = pTerm * termFactor
            elif termFactor.is_Pow or termFactor.is_Symbol:
                # Handling of x^y and x^1.
                base_exp = termFactor.as_base_exp() if termFactor.is_Pow else (termFactor, 1)
                # Retrieve the indices of the coefficient.
                ix = base_exp[0].name[base_exp[0].name.find('_'):]
                n = int(ix.split('_')[1])
                k = int(ix.split('_')[2])
                # Build the a_nk sum.
                pTerm = pTerm * a_nk(n, k) ** base_exp[1]
            else:
                raise Exception("Unhandled factor in sympy expression:" + str(termFactor))
        pTot = pTot + pTerm
    return pTot
