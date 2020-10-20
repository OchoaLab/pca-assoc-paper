###################
### SIMULATIONS ###
###################

# main simulations follow the real dataset analysis in parallel, so that we take advantage of its embarrasingly parallel structure (for potential cluster runs)

# all have: -k 10 -f 0.1
# this is like my older FST work
# but unlike my main gas-rgls sim (instead has: -k 3 -f 0.3)

### construct admixture proportions, ancestral inbreeding, and family structure if needed
# done once, shared across replicates

# the large sample size simulation
time Rscript sim-00-sim-pop.R # -n 1000 -g 1
# 0m1.047s ideapad

# the small sample size simulation
time Rscript sim-00-sim-pop.R -n 100 # -g 1
# 0m0.783s ideapad

# the family structure simulation
time Rscript sim-00-sim-pop.R -g 20 # -n 1000
# 0m13.656s ideapad (concurrent with plink, so probably would have been faster)

### construct genotype matrices

# params shared across reps
# use right set for each case

## LARGE SAMPLE SIZE
# sample size (normal n=1000, small n=100)
n=1000
# number of causal loci should scale with n
mc=$(( n / 10 ))
# to create family structure (g=20) or not (g=1)
g=1
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"

## FAMILY STRUCTURE
# sample size (normal n=1000, small n=100)
n=1000
# number of causal loci should scale with n
mc=$(( n / 10 ))
# to create family structure (g=20) or not (g=1)
g=20
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"

## SMALL SAMPLE SIZE
# sample size (normal n=1000, small n=100)
n=100
# number of causal loci should scale with n
mc=$(( n / 10 ))
# to create family structure (g=20) or not (g=1)
g=1
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"

for rep in {1..50}; do
    # draw random genotypes
    time Rscript sim-01-draw-geno.R -r $rep -n $n -g $g
    # 0m19.112s ideapad (concurrent with plink)

    # draws a random trait
    # only place where --m_causal gets specified (not in file paths, etc)
    time Rscript sim-02-sim-trait.R --bfile $name -r $rep --m_causal $mc
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

    # PCA runs (with plink)
    # same but with standard PCA estimates (from my R code, kinship_std ROM version)
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --std
    # 0m3.422s ideapad (concurrent with plink)
    # do all PCs of this rep
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs
    done
    # cleanup
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean --std
    
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

# creates final plot for paper!
time Rscript real-09-figs.R --bfile $name
# exploratory version compares pure PCA to PCs from popkinsuppl::kinship_std (practically the same)
time Rscript real-09-figs.R --bfile $name --pca

# a summary of "best PCs" in a simple analysis
Rscript real-13-stats.R --bfile $name
# LARGE
#   method         metric  best   min
# 1 pca-plink-pure rmsd      45     3
# 2 pca-plink-pure auc        3     3
# 3 gcta           rmsd       0     0
# 4 gcta           auc        1     0
# best rmsd: pca-plink-pure (significant)
# best auc: gcta (significant)
#
# FAMILY
#   method         metric  best   min
# 1 pca-plink-pure rmsd      90    84
# 2 pca-plink-pure auc       10     4
# 3 gcta           rmsd       0     0
# 4 gcta           auc        0     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)
#
# SMALL
#   method         metric  best   min
# 1 pca-plink-pure rmsd      88     3
# 2 pca-plink-pure auc        2     1
# 3 gcta           rmsd       0     0
# 4 gcta           auc        0     0
# best rmsd: pca-plink-pure (tie)
# best auc: gcta (tie)


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

# a comparison of RMSD and lambda across all datasets
time Rscript real-11-inflation-across-datasets.R

# archive p-values and individual summary files (move out of space that gets synced between computers; these files take up loads of space and are no longer needed)
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90
