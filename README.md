# Point-to-ellipse series

Python (reference) and C++ (fast) implementations of the series expansions and
rational coefficients for the point-to-ellipse relation:

**Part I — expansion far from the center:**
Nilsson, J.-O. *Point-to-ellipse Fourier series.*
arXiv: [2507.08807](https://doi.org/10.48550/arXiv.2507.08807) \[math.GM]
(a copy is included: [2507.08807v1.pdf](2507.08807v1.pdf)).

**Part II — expansion about the center (inside-evolute series):** *in preparation.*

## Introduction

The point-to-ellipse relation — the normal direction φ and the signed distance
(height) h from an ellipse to a nearby point — is normally computed by iteration
or by solving the quartic directly. These works instead give general **series
expansions**: Part I valid far from the center, Part II inside the ellipse's
evolute. This repository provides their implementations: exact symbolic
computation of the rational coefficients and evaluation of the truncated series.
For the mathematics, see the articles.

## Python (reference implementation)

`python/main.py` demonstrates the series by converting a reference
latitude/altitude to Cartesian coordinates and reconstructing the
point-to-ellipse relations from the truncated series (multi-precision).

- `expansions.py`, `expansions_evo.py` — the truncated series expansions
  (far-field / inside-evolute; sin-power form, with Fourier multiple-angle variants).
- `coefficients.py`, `coefficients_evo.py` — the exact rational series
  coefficients (far-field / inside-evolute).
- `polynomials.py`, `series_substitutions.py`, `series.py` — the symbolic
  machinery (potential/Bell polynomials, series arithmetic, substitutions).
- `ellipse.py` — the closed-form (Vermeille) reference and ellipse parameters.
- `symbols.py`, `cache.py`, `util.py` — symbols, memoization, utilities.
- `analysis/` — convergence plots; `tests/` — pytest checks against the
  closed-form reference.

The Python code is a readable reference (practical to ~order 15).

## C++ (fast implementation)

`cpp/` mirrors the Python coefficients and series — checked against the Python
reference via CSV in the tests — but reaches far higher orders. Coefficients are
exact GMP rationals; the series are templated and evaluate either symbolically
(SymEngine) or numerically (`mpfr::mpreal`).

- `coefficients[_evo].{hpp,cpp}`, `expansions[_evo].hpp` — coefficients and series.
- `convergence/` — `sweep_evo`, `diag_coeff_evo`: empirical convergence / region-of-
  convergence studies (CSV output, plotted by `python/analysis`).
- `benchmarks/`, `tests/` — micro-benchmarks and Python-cross-checked tests.

With the C++ code, the practical order is >40.

Build (CMake ≥ 3.15, C++20):

    cmake -S cpp -B cpp/build && cmake --build cpp/build

## Dependencies

- **Python:** sympy, mpmath (matplotlib, numpy for the analysis plots).
- **C++:** a C++20 compiler, CMake ≥ 3.15, GMP/gmpxx, MPFR, SymEngine.

## License

BSD 2-Clause (see [LICENSE.txt](LICENSE.txt)).

## Contact

john_nil at hotmail period com.
