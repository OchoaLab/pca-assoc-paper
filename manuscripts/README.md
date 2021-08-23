# PCA association manuscript

This folder contains the manuscript in PDF format and the various LaTeX source files.
However, figures and tables required for this manuscript are available in `../data/`.

The LaTeX source depends on many standard packages available on texlive and other distributions.
The only exception is [kinshipsymbols](https://github.com/ochoalab/kinshipsymbols).

To compile manuscript from scratch, you need to run these commands (several times as shown since references need to be updated):
```bash
pdflatex pca-assoc
bibtex pca-assoc
pdflatex pca-assoc
pdflatex pca-assoc
```
