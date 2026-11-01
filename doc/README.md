# Manuscript sources

The root-level PDFs are versioned snapshots: `2507.08807v1.pdf` is Part I's
arXiv v1, and `point2ellipse_inverse_v0.pdf` is Part II's submitted version,
pending replacement with the arXiv PDF once available.

The sources in this directory may contain subsequent corrections and need not
match those snapshots. Part I already includes minor corrections, with a list
of edits in its source comments. Corrections may accumulate before a revised
version is submitted to arXiv.

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
