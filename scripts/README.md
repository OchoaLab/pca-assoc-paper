## PCA association project code

### Dependencies

To install R package dependencies see `packages.R`.

Binary dependencies are [plink2](https://www.cog-genomics.org/plink/2.0/), [GCTA](https://yanglab.westlake.edu.cn/software/gcta/), and [Eigensoft](https://github.com/DReichLab/EIG).

### Reproducing results

The [real data processing instructions](https://github.com/OchoaLab/data) are in the "data" repository because shared with other projects.
Those contain complete instructions for obtaining public datasets and producing the processed versions that are the inputs to our analysis.

The pipeline for reproducing our results is described in the `main.bash` file.

However, the evaluations are very computationally intensive, particularly the ones involving the real datasets, which cannot be run on a desktop in a reasonable amount of time.
We used computer clusters to perform most of the plink2 and GCTA association runs.
Also, there are ocasional use of absolute paths that only make sense in our specific computers.
Therefore, discretion is advised as to how to best adapt and run this code in your given setup.
