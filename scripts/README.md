## PCA association project code

The pipeline for reproducing our results is described in the two `main-*.bash` files.
The `main-sim.bash` commands were run first, followed by `main-real.bash`.

However, the evaluations are very computationally intensive, particularly the ones involving the real datasets, which cannot be run on a desktop within a reasonably short amount of time.
We used computer clusters to perform most of the plink2 and GCTA association runs.

