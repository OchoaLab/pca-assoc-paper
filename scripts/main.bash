# shared parameters

# non-default low-heritability simulation
h=0.3
mcf=27 # 10*8/3 rounded, adjusts expected coefficient size for decrease in heritability
# and env variances (only used with low herit)
env1=0.3
env2=0.2


###################
### SIMULATIONS ###
###################

## LARGE SAMPLE SIZE
# sample size (normal n=1000, small n=100)
n=1000
# to create family structure (g=20) or not (g=1)
g=1
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"
# run shared steps
. main-sim.bash


## SMALL SAMPLE SIZE
# sample size (normal n=1000, small n=100)
n=100
# to create family structure (g=20) or not (g=1)
g=1
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"
# run shared steps
. main-sim.bash


## FAMILY STRUCTURE
# sample size (normal n=1000, small n=100)
n=1000
# to create family structure (g=20) or not (g=1)
g=20
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"
# run shared steps
. main-sim.bash


######################
### REAL GENOTYPES ###
######################

# creates a table other scripts might use, especially when creating figures, for consistent human-facing naming, and loops across all datasets
# writes data/datasets.txt
Rscript real-00-datasets.R


# shared by R script calls further below
name='HoPacAll_ld_prune_1000kb_0.3_maf-0.01'
DATA_DIR='/home/viiia/dbs/humanOrigins'
# run shared steps (real, tree simulations, and 4th degree relatives removed)
. main-real.bash


# version for HGDP WGS
name='hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1'
DATA_DIR='/home/viiia/dbs/hgdp_wgs'
# run shared steps
. main-real.bash


# version for TGP
name='tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01'
DATA_DIR='/home/viiia/dbs/tgp-nygc'
# run shared steps
. main-real.bash


###############
### GLOBALS ###
###############

# cross-dataset analyses, including all of: sim, real, and real-sim
# this sequence also creates all of the main figures/tables in order of appearance

# reports actual dimensions (n_ind, m_loci, K, and m_causal) used in all datasets
# (validates every replicate too! for m_causal compares both FES and RC!)
# writes data/dimensions.txt
time Rscript real-14-dimensions.R
# 0m17.481s viiiaR5

# MAF plot code
time Rscript real-17-mafs-plot.R
# 0m13.937s viiiaR5

# popkin kinship plots
# includes data generated by MAF script above
time Rscript all-01-kinship.R
# 30s ideapad

# a simple figure that illustrates the methods
# (uses data from "Admix. Large sim." only)
time Rscript sim-10-measures-fig.R 2
# Inflation factor: pca-plink-pure: 2.62753790231827
# Inflation factor: gcta: 0.890044438627719
# Inflation factor: gcta: 1.00079570634139

# main statistical evaluations between methods
time Rscript real-13-stats.R
# Number of tests (for Bonferroni): 72
time Rscript real-13-stats.R --herit $h --m_causal_fac $mcf
# Number of tests (for Bonferroni): 48
time Rscript real-13-stats.R --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l
# Number of tests (for Bonferroni): 72

# final plot gathers all three datasets into a single multipanel figure
time Rscript real-15-plots-big.R
time Rscript real-15-plots-big.R --fes
time Rscript real-15-plots-big.R --real
time Rscript real-15-plots-big.R --real --fes
time Rscript real-15-plots-big.R --real_sim
time Rscript real-15-plots-big.R --real_sim --fes
# low herit versions
time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf
time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --fes
time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --real
time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --real --fes
# these two don't exist
# time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --real_sim
# time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --real_sim --fes
# env versions (all with -l to include gcta-labs version)
time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l
time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l --fes
time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l --real
time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l --real --fes
# these two don't exist
# time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l --real_sim
# time Rscript real-15-plots-big.R --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l --real_sim --fes


# calculages popkin and eigensoft eigenvectors, calculates TW stats, makes plot
time Rscript all-02-eigen.R
# 30m14.103s viiiaR5 first time

# calculate numbers of significant PCs against trait (null model without SNPs)
# only ran these two versions, to match the data that was there when "eigen" is presented (i.e. low herit and env come later, but eigen stuff doesn't come back at that point; don't expect major differences)
# (also, will have to exclude real-sims for those cases, current code doesn't and will just die)
time Rscript all-05-pcs-num-sig.R
# 5m4.593s
time Rscript all-05-pcs-num-sig.R --fes

# a comparison of RMSD and lambda across ALL datasets (including FES and RC traits)
# this version that fits top half only (makes most sense for our goal of talking mostly about inflation)
# same script now also compares RMSD and type I error, and also AUC vs power
time Rscript real-11-inflation-across-datasets.R
# model fit:
# rmsd ~ a * (lambda^b - 1) / (lambda^b + 1)
#         a         b 
# 0.5644519 0.6187852 
# threshold map (sigmoidal): lambda = 1.05, RMSD = 0.00852
# Inverse threshold map (sigmoidal): RMSD = 0.01, lambda = 1.06

# create figure that compares h3 and env simulations directly
time Rscript all-04-plot-low-herit-vs-env.R

# look at kinship distributions in real and simulated datasets
time Rscript all-03-king.R

# calculate dimensions table for this subset separately
time Rscript king-00-dimensions.R

# make AUC/RMSD plot and stats tests for the limited king-cutoff test
time Rscript king-01-rmsd-auc-plot.R
# HoPacAll_ld_prune_1000kb_0.3_maf-0.01_king-cutoff-4
# RMSD: LMM r = 0
# AUC: LMM r = 0
# hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1_king-cutoff-4
# RMSD: LMM r = 0
# AUC: LMM r = 0
# tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01_king-cutoff-4
# RMSD: LMM r = 0
# AUC: LMM r = 0

# same for low herit tests
time Rscript king-01-rmsd-auc-plot.R --herit $h --m_causal_fac $mcf
# HoPacAll_ld_prune_1000kb_0.3_maf-0.01_king-cutoff-4
# RMSD: LMM r = 0
# AUC: LMM r = 0
# hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1_king-cutoff-4
# RMSD: LMM r = 0
# AUC: LMM r = 0
# tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01_king-cutoff-4
# RMSD: LMM r = 0
# AUC: LMM r = 0

# and env, here add gcta-labs!
time Rscript king-01-rmsd-auc-plot.R --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l
# HoPacAll_ld_prune_1000kb_0.3_maf-0.01_king-cutoff-4
# RMSD: LMM r = 10
# AUC: LMM lab. r = 0
# hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1_king-cutoff-4
# RMSD: LMM r = 0/PCA r = 20
# AUC: LMM lab. r = 0
# tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01_king-cutoff-4
# RMSD: LMM lab. r = 0
# AUC: LMM lab. r = 0

