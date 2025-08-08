
""" Various Utility functions for point-to-ellipse series expansions.
"""

# External includes.
from collections.abc import Callable
import sympy as sp
from functools import lru_cache

# Internal includes.
import polynomials
from symbols import e2
import series


@lru_cache(maxsize=None)
def poly_bell_substitution(p : sp.core.Expr, b : str) -> series.SeriesBase:
    """Expand a polynomial p(a_0,...,a_n) by substituting Bell polynomials for a_n^i.

    p is a multidimensional polynomial in variables a_n. (The assumption is
    that the symbol name ends with "_<number>" which can be parsed and
    transferred to the substitution variable.)
    Assume

    a_n = \sum_{k=1}^\infty b_{n,k} z^k

    where b_{n,k} are new variables, then

    a_n^i = \sum_{k=j} \hat{B}_{k,i}(b_{n,1},...,b_{n,k-i+1})z^k

    where \hat{B}_{k,i}(...) is a partial ordinary Bell polynomial.
    See https://en.wikipedia.org/wiki/Bell_polynomials

    This function makes the Bell polynomial substitution.

    :param p: The polynomial to make the substitution in. Polynomial must be expanded.
    :param b: The variable name of series coefficients.
    :return: A coefficient sequence representing the substituted polynomial.
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
                # Retrieve the indices of the coefficient.
                ix = base_exp[0].name[base_exp[0].name.find('_'):]
                def bellLambdaGen(j, x):
                    return lambda n : polynomials.partial_ordinary_bell_polynomial(n, j, x) if n >= j else 0
                seqTerm = seqTerm * series.Series(bellLambdaGen(base_exp[1], b + ix))
            else:
                raise Exception("Unhandled factor in sympy expression:" + str(termFactor))
        seqTot = seqTot + seqTerm
    return seqTot

@lru_cache(maxsize=None)
def a_nk_sub(p : sp.core.Expr, n_offset : int, d_nkl : Callable[[int, int, int], sp.core.Expr]) -> sp.core.Expr:
    """Substitute a_{n,k} for \sum_{l=\max(k,n+n_offset)}^{n+k} d_{n,k,l} e²^l in a polynomial.

    a_{n,k} is the arguments of the Bell polynomial.

    :param d_nkl: The coefficient of the substitution sum.
    :param n_offset: Offset in n in the lower limit of the sum.
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
                a_nk = sp.Integer(0)
                for l in range(max(k, n + n_offset), n + k + 1):
                    a_nk = a_nk + d_nkl(n, k, l) * e2 ** l
                pTerm = pTerm * a_nk ** base_exp[1]
            else:
                raise Exception("Unhandled factor in sympy expression:" + str(termFactor))
        pTot = pTot + pTerm
    return pTot

def sin_pow_to_cos_mul(n : int, k : int, l : int, d_nkl : Callable[[int, int, int], sp.core.Expr],
                       n_min : int, k_pp : int) -> sp.core.Rational | sp.core.Integer:
    """Fourier multiple-angle cos series coefficient from sin-power series.

    :param n: sin-multiple.
    :param k: rho power.
    :param l: e³ power.
    :param d_nkl: Sin-power series coefficients.
    :param n_min: Lowest sin-power.
    :param k_pp: Maximum e² power offset in sin-power series.
    :return: Fourier multiple-angle cos series coefficient of
    """
    c_nkl = sp.Integer(0)
    for i in range(max(max(n, n_min), l - k - k_pp), l + 1):
        tmp = d_nkl(i, k, l) * sp.Rational(2, 2 ** (2 * i)) * (-1) ** n * sp.binomial(2 * i, i - n)
        c_nkl = c_nkl + (tmp if n != 0 else tmp / 2)
    return c_nkl
