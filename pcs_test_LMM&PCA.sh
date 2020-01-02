#!/bin/bash
#
#SBATCH --job-name=pcs_test_LMM&PCA
#SBATCH --mail-user=yy222@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=8000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=pcs_test_LMM&PCA.stdout
#SBATCH --error=pcs_test_LMM&PCA.stderr


Rscript pcs_test_LMM&PCA.R

