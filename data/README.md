## PCA association project data

This folder holds some raw data (raw tables with SRMSD and AUCPR values per dataset, number of PCs, and replicates) as well as other summaries, including main tables and figures used in the manuscript.

There are 9 datasets, each in a subfolder with their own dataset-specific data.
The mapping between directory names and the names used in the paper is given in `datasets.txt`.
Each dataset has results for the "random coefficients" RC trait model in their base, and corresponding results for the "fixed effect sizes" (FES) trait models in their respective `fes/` subdirectories.
Likewise, main figures for RC traits are in the current directory, and for FES traits they are in the `fes/` subdirectory.

Raw genotype data is not available in this repository, due to its size.
The real datasets (1000 Genomes, HGDP, and Human Origins) are available in full elsewhere, see [real data processing instructions](https://github.com/OchoaLab/data) for download and pre-processing instructions.
All of the genotype and phenotype simulations were performed with publicly-available software, including data and code provided in this repository.
