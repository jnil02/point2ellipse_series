
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
    """Expand a polynomial p(a_0,...,a_n) by substituting Bell polynomials for a_n^i.

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
    """Coefficient of the power of a double power series where the first series start from 0 and the second starts from 1.

    :param n: First index of resulting series coefficients.
    :param i: The power of the double power series.
    :return: Series representing the coefficient.
    """
    # Polynomial for b_{n,i} in terms of {a_0,...,a_n}.
    b_ni = polynomials.ordinary_potential_polynomial(n, i, "a")
    # Polynomial for the varrho^k coefficient in b_{n,i} in terms of {a_{n,1},...a_{n,k+1}}
    return poly_bell_substitution(b_ni)

@cache.ints_cache
def double_series_power_coeff_evo(n: int, i: int) -> series.SeriesBase:
    """Coefficient of the power of a double power series for the inside-evolute chain.

    Same as double_series_power_coeff, but the inner generator a_l starts at
    z^{l+1} rather than z^1. This tightens the Bell substitution's nonzero
    condition to k>=i*(l+1) (manuscript: the kappa-sum defining c^phi_{k,l,n,i}
    starts at i_l(l+1)). The result is unchanged from the untightened version;
    the pruned terms are identically zero after the a_{n,k} substitution.

    :param n: First index of resulting series coefficients.
    :param i: The power of the double power series.
    :return: Series representing the coefficient.
    """
    b_ni = polynomials.ordinary_potential_polynomial(n, i, "a")
    return poly_bell_substitution(b_ni, start=lambda l: l + 1)

def a_nk_ser(n: int, k: int, n_offset: int, d_nkl: Callable[[int, int, int], sp.core.Expr], e2: sp.core.Symbol) -> sp.core.Expr:
    """Specific finite a_{n,k} series from max(k, n+offsetI to n+k.

    :param n: n index
    :param k: k index
    :param n_offset: offset
    :param d_nkl: tripple sum coefficient.
    :param e2 series base variable.
    :return: expression for finite series.
    """
    a_nk = sp.S.Zero
    for l in range(max(k, n + n_offset), n + k + 1):
        a_nk = a_nk + d_nkl(n, k, l) * e2 ** l
    return a_nk

def a_nk_C(n: int, k: int, c: Callable[[int, int, int], sp.core.Expr], e2: sp.core.Symbol):
    """Specific finite a_{n,k} series from 1 to k if k>=n+1.

    :param n: n index
    :param k: k index
    :param c: tripple sum coefficient.
    :param e2 series base variable.
    :return: expression for finite series.
    """
    a_nk = sp.S.Zero
    if k < n+1:
        return a_nk
    for l in range(1, k + 1):
        a_nk += c(n, k, l) * e2 ** l
    return a_nk

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
