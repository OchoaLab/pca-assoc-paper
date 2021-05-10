###################
### SIMULATIONS ###
###################

# main simulations follow the real dataset analysis in parallel, so that we take advantage of its embarrasingly parallel structure (for potential cluster runs)

# all have: -k 10 -f 0.1
# this is like my older FST work
# but unlike my main gas-rgls sim (instead has: -k 3 -f 0.3)

# params shared across reps
# use right set for each case

## LARGE SAMPLE SIZE
# sample size (normal n=1000, small n=100)
n=1000
# to create family structure (g=20) or not (g=1)
g=1
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"

## SMALL SAMPLE SIZE
# sample size (normal n=1000, small n=100)
n=100
# to create family structure (g=20) or not (g=1)
g=1
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"

## FAMILY STRUCTURE
# sample size (normal n=1000, small n=100)
n=1000
# to create family structure (g=20) or not (g=1)
g=20
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"

# construct admixture proportions, ancestral inbreeding, and family structure if needed
# done once, shared across replicates
time Rscript sim-00-sim-pop.R -n $n -g $g
# large: 0m1.047s ideapad
# small: 0m0.783s ideapad
# family: 0m13.656s ideapad (concurrent with plink, so probably would have been faster)

### construct genotype matrices, traits, etc

for rep in {1..50}; do
    # draw random genotypes
    time Rscript sim-01-draw-geno.R -r $rep -n $n -g $g
    # 0m19.112s ideapad (concurrent with plink)

    # draws a random trait
    time Rscript sim-02-sim-trait.R --bfile $name -r $rep
    # 0m1.914s ideapad (concurrent with plink)

    # preprocess with GCTA (makes GRM and max PCs)
    time Rscript real-00-preprocess-gcta.R --bfile $name/rep-$rep
    # 0m56.426s ideapad (concurrent with plink)

    # get PCs in R, using my formula
    # NOTE: uses ROM version (known bias, does not match popular versions; I think this is best)
    # top 4 PCs match perfectly here, 5 agrees partially, after that it's a mess
    time Rscript real-01-pca-test.R --bfile $name/rep-$rep
    # 0m17.695s ideapad (concurrent with plink)

    # get PCs using plink2
    time Rscript real-01-pcs-plink.R --bfile $name/rep-$rep
    
    # NOTE: unlike real dataset, here we don't need popkin estimates, we use true p_anc instead (for simulating trait earlier)!

    # GCTA runs
    # this creates auxiliary GCTA PCA files (redundant, will delete when done with this analysis)
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep
    # 0m3.215s ideapad (concurrent with plink)
    # do all PCs of this rep
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs
    done
    # cleanup
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean

    # # PCA runs (with plink)
    # # same but with standard PCA estimates (from my R code, kinship_std ROM version)
    # time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --std
    # # 0m3.422s ideapad (concurrent with plink)
    # # do all PCs of this rep
    # for pcs in {0..90}; do
    # 	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs
    # done
    # # cleanup
    # time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean --std
    
    # PCA runs (with *pure* plink)
    # lastly, same but with PCs from plink2
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink
    # do all PCs of this rep
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs --plink
    done
    # cleanup
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink --clean
done

# for MAF comparisons plot (only needs rep-1)
time Rscript real-16-mafs.R --bfile $name --sim


# summarizes p-values into AUC and RMSD for each method/rep/pc
time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90
# 46m10.111s viiiaX6 (large)
# 38m9.266s viiiaX6 (family)
# 33m3.880s viiiaX6 (small)
# 11m37.970s ideapad (small plink pure)
# 12m18.912s ideapad (large plink pure)
# 9m58.027s + 2m49.465s ideapad (family plink pure)

# read all individual summary tables (tiny files), gather into a master table
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
# 0m32.340s viiiaX6 (large)
# 0m32.287s viiiaX6 (family)
# 0m32.228s viiiaX6 (small)
# 0m38.671s ideapad (large plink pure)
# 0m36.665s ideapad (family plink pure)
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 -a # to run on fully-archived data

# creates final plot for paper!
time Rscript real-09-figs.R --bfile $name
# exploratory version compares pure PCA to PCs from popkinsuppl::kinship_std (practically the same)
time Rscript real-09-figs.R --bfile $name --pca

# a summary of "best PCs" in a simple analysis
Rscript real-13-stats.R --bfile $name
# LARGE
#   method         metric  best   min
# 1 pca-plink-pure rmsd      21     3 # OLD 45,3
# 2 pca-plink-pure auc        3     3
# 3 gcta           rmsd       0     0
# 4 gcta           auc        1     0
# best rmsd: pca-plink-pure (tie) # OLD significant
# best auc: gcta (significant)
#
# SMALL
#   method         metric  best   min
# 1 pca-plink-pure rmsd       5     2 # OLD 88,3
# 2 pca-plink-pure auc        2     1
# 3 gcta           rmsd       0     0
# 4 gcta           auc        0     0
# best rmsd: gcta (significant) # OLD: pca-plink-pure (tie)
# best auc: gcta (tie)
#
# FAMILY
#   method         metric  best   min
# 1 pca-plink-pure rmsd      90    84
# 2 pca-plink-pure auc       10     4
# 3 gcta           rmsd       0     0
# 4 gcta           auc        0     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)

# tests that p-value vectors have the right lengths of m_loci
# to make sure nothing was corrupted due to scripts stopping unexpectedly or incomplete file transfers
# (now test is peformed within real-07-auc-rmsd.R above too, but this can retest everything without recalculating expensive summary statistics).
# --final requires that all files exist!
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --final
# 2m54.797s - 3m43.697s viiiaX6
# 4m29.032s ideapad (small plink pure)
# 2m32.743s ideapad (large plink pure)
# 4m51.305s ideapad (family plink pure)
###

# archive p-values and individual summary files (move out of space that gets synced between computers; these files take up loads of space and are no longer needed)
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90

# a simple figure that illustrates the methods
Rscript sim-10-measures-fig.R 13
# Inflation factor: gcta: 0.977811992816224
# Inflation factor: pca-plink-pure: 2.55503828726139
# Inflation factor: gcta: 0.861113356022774

# final plot gathers all three simulations into a single multipanel figure
time Rscript real-15-plots-big.R

########################
### const_herit_loci ###
########################

# addition to complement existing analysis with alternative trait from const_herit_loci model
# not sure if it will make a meaningful difference
# NOTE: many steps that depend on genotypes only aren't redone (are shared from prev run)

for rep in {1..50}; do
    time Rscript sim-02-sim-trait.R --bfile $name -r $rep --const_herit_loci

    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs --plink --const_herit_loci
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink --clean

    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs --const_herit_loci
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean
done

time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90 --const_herit_loci
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --const_herit_loci

time Rscript real-09-figs.R --bfile $name --const_herit_loci
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --final --const_herit_loci
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --const_herit_loci -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --const_herit_loci

Rscript real-13-stats.R --bfile $name --const_herit_loci
# LARGE
#   method         metric  best   min
# 1 pca-plink-pure rmsd      90     8 # OLD: 90,4
# 2 pca-plink-pure auc        4     3
# 3 gcta           rmsd       0     0
# 4 gcta           auc        1     0
# best rmsd: gcta (tie) # OLD: pca-plink-pure (significant)
# best auc: gcta (significant)
#
# SMALL
#   method         metric  best   min
# 1 pca-plink-pure rmsd       3     2 # OLD 84,2
# 2 pca-plink-pure auc        1     1
# 3 gcta           rmsd       0     0
# 4 gcta           auc        0     0
# best rmsd: gcta (significant) # OLD: tie
# best auc: gcta (significant)
#
# FAMILY
#   method         metric  best   min
# 1 pca-plink-pure rmsd      90    85
# 2 pca-plink-pure auc       24     6 # OLD: 24,5
# 3 gcta           rmsd       0     0
# 4 gcta           auc        0     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)

time Rscript real-15-plots-big.R --const_herit_loci

