# Manuscript sources

The root-level PDFs are versioned snapshots:
[2507.08807v1.pdf](../2507.08807v1.pdf) is Part I's arXiv v1, and
[2610.40690v1.pdf](../2610.40690v1.pdf) is Part II's arXiv v1.

The sources in this directory may contain subsequent corrections and need not
match those snapshots. Part I already includes minor corrections, with a list
of edits in its source comments. Corrections may accumulate before a revised
version is submitted to arXiv.

## Series implementations

Each manuscript presents its series expansions in a single theorem. The
corresponding truncated series and exact rational coefficients are implemented
in the following files:

| Manuscript | Python series / coefficients | C++ series / coefficients |
| --- | --- | --- |
| Part I: expansion far from the center | [expansions.py](../python/expansions.py), [coefficients.py](../python/coefficients.py) | [expansions.hpp](../cpp/expansions.hpp), [coefficients.cpp](../cpp/coefficients.cpp) |
| Part II: expansion about the center | [expansions_evo.py](../python/expansions_evo.py), [coefficients_evo.py](../python/coefficients_evo.py) | [expansions_evo.hpp](../cpp/expansions_evo.hpp), [coefficients_evo.cpp](../cpp/coefficients_evo.cpp) |

The numerical tests check these implementations; they are not numerical results
reported in the manuscripts. See the manuscripts for the formulas and their
mathematical conventions.

## Building

Install a LaTeX distribution providing pdfLaTeX, Biber, and the packages used
by the sources, including `amsart`, `biblatex`, TikZ, `caption`, `mathtools`,
`tabularx`, `ragged2e`, `hyperref`, and `longtable`.

Run each set of commands from the repository root. The working directory must
be the manuscript's folder so that bibliography, figure, and table paths resolve.

Part I:

```bash
cd doc/part_I_expansion_about_inf
pdflatex -interaction=nonstopmode -halt-on-error point2ellipse_series.tex
biber point2ellipse_series
pdflatex -interaction=nonstopmode -halt-on-error point2ellipse_series.tex
pdflatex -interaction=nonstopmode -halt-on-error point2ellipse_series.tex
```

Part II (start again from the repository root):

```bash
cd doc/part_II_expansion_about_zero
pdflatex -interaction=nonstopmode -halt-on-error point2ellipse_inverse.tex
biber point2ellipse_inverse
pdflatex -interaction=nonstopmode -halt-on-error point2ellipse_inverse.tex
pdflatex -interaction=nonstopmode -halt-on-error point2ellipse_inverse.tex
```

These commands produce a PDF and auxiliary files in the corresponding source
folder. They do not update the root-level PDF snapshots. The figure PDFs and
coefficient tables needed by both manuscripts are already included.
