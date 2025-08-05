
# Point-to-ellipse Fourier series

This repository contains a naive Python implementation of series and
coefficients of:

<b>Nilsson, J.-O. Point-to-ellipse Fourier series. doi:
[https://doi.org/10.48550/arXiv.2507.08807](https://doi.org/10.48550/arXiv.2507.08807)
arXiv: 2507.08807 [math.ag]</b>

<i>Note, a copy of the article is found here in the repo,
[2507.08807v1.pdf](2507.08807v1.pdf).</i>

The code may be used to compute coefficients up to order ~15. Beyond that it
simply takes too long time. A refined (C++) implementation could probably
provide significantly more coefficients.

## Introduction

Point-to-ellipse relations are normally computed with iterative methods or by
directly solving the quartic latitude equation. In the article, the first
general series expansions for the point-to-ellipse relations are provided. The
series expansions are interesting for analysis and provide many computational
possibilities. Here their actual implementations are provided, i.e. symbolic
computations of coefficients and series up to some specified order. For
mathematical details, see the article.

## The code

The main entry-point is the main.py which show-case the series by computing a
number of point-to-ellipse relations. The actual implementations are found in
the following packages:

- series.py -- Symbolic implementations of the (truncated) Fourier and sin-power series.
- coefficients.py -- Computations of rational series coefficients.
- ellipse.py -- Ellipse parameters for which the series are computed and some related utility functions.
- util.py -- Various utility functions.
- symbols.py -- Definition of sympy symbols used in the code.

## Dependencies

The following external packages are imported
 - sympy
 - mpmath
 - functools

## License

BSD 2-Clause License

## Contact
john_nil at hotmail period se.