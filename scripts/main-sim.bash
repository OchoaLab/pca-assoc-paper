
# construct admixture proportions, ancestral inbreeding, and family structure if needed
# done once, shared across replicates
time Rscript sim-00-sim-pop.R -n $n -g $g
# large: 0m1.047s ideapad
# small: 0m0.783s ideapad
# family: 0m0.440s viiiaR5

### construct genotype matrices, traits, etc

for rep in {1..50}; do
    # draw random genotypes
    time Rscript sim-01-draw-geno.R -r $rep -n $n -g $g
    # 0m19.112s ideapad (concurrent with plink)
    # 0m50.744s viiiaR5 family 26% max mem with Rcpp!

    # draws a random trait
    time Rscript real-04-simtrait.R --bfile $name -r $rep --sim

    # preprocess with GCTA (makes GRM and max PCs)
    time Rscript real-00-preprocess-gcta.R --bfile $name/rep-$rep

    # get PCs using plink2
    time Rscript real-01-pcs-plink.R --bfile $name/rep-$rep
    
    # NOTE: unlike real dataset, here we don't need popkin estimates, we use true p_anc instead (for simulating trait earlier)!

    # GCTA runs
    # this creates auxiliary GCTA PCA files (redundant, will delete when done with this analysis)
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep
    # do all PCs of this rep
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs
    done
    # cleanup
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean

    # PCA runs (with *pure* plink)
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink
    # do all PCs of this rep
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs
    done
    # cleanup
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink --clean
done

# estimate local kinship with KING-robust, rep-1 only
cd ../data/$name/rep-1
time plink2 --bfile data --make-king triangle bin4 --out data
# cleanup
rm data.log 
# return to scripts dir
cd ../../../scripts/

# for MAF comparisons plot (only needs rep-1)
time Rscript real-16-mafs.R --bfile $name --sim
# get popkin estimates for an overview plot only (also rep-1/ only)
time Rscript real-03-popkin.R --bfile $name --sim


# summarizes p-values into AUC and RMSD for each method/rep/pc
time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90
# 11m37.970s ideapad (small)
# 12m18.912s ideapad (large)
# 4m32.344s viiiaR5 family (first run parallelized)

# read all individual summary tables (tiny files), gather into a master table
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
# 0m38.671s ideapad (large)
# 0m36.665s ideapad (family)
# 0m48.798s viiiaR5 family (slower old HD?)
###time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 -a # to run on fully-archived data

# creates final plot for paper!
time Rscript real-09-figs.R --bfile $name

# tests that p-value vectors have the right lengths of m_loci
# to make sure nothing was corrupted due to scripts stopping unexpectedly or incomplete file transfers
# (now test is peformed within real-07-auc-rmsd.R above too, but this can retest everything without recalculating expensive summary statistics).
# --final requires that all files exist!
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --final
# 4m29.032s ideapad (small)
# 2m32.743s ideapad (large)
# 1m5.998s viiiaR5 family

# archive p-values and individual summary files (move out of space that gets synced between computers; these files take up loads of space and are no longer needed)
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90

###########
### FES ###
###########

# addition to complement existing analysis with trait from FES model
# NOTE: steps that depend only on genotypes aren't redone (are shared from prev run)

for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep --sim --fes

    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs --fes
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink --clean

    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs --fes
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean
done

time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90 --fes
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --fes
time Rscript real-09-figs.R --bfile $name --fes
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --final --fes
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --fes -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --fes

#################
### herit 0.3 ###
#################

h=0.3
mcf=27 # 10*8/3 rounded, adjusts expected coefficient size for decrease in heritability

# do both RC and FES, instructions together here
for rep in {1..50}; do
    # trait simulations were performed locally
    time Rscript real-04-simtrait.R --bfile $name -r $rep --sim --herit $h --m_causal_fac $mcf
    time Rscript real-04-simtrait.R --bfile $name -r $rep --sim --herit $h --m_causal_fac $mcf --fes

    # PCA and LMM runs were actually performed on DCC (see batch.R)
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs --herit $h --m_causal_fac $mcf
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs --herit $h --m_causal_fac $mcf --fes
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink --clean

    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs --herit $h --m_causal_fac $mcf
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs --herit $h --m_causal_fac $mcf --fes
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean
done

# RC
time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf
time Rscript real-09-figs.R --bfile $name --herit $h --m_causal_fac $mcf
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --final
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf

# FES
time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes
time Rscript real-09-figs.R --bfile $name --herit $h --m_causal_fac $mcf --fes
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --final --fes
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes

###########
### ENV ###
###########

# start from prev low-herit sim
h=0.3
mcf=27 # 10*8/3 rounded, adjusts expected coefficient size for decrease in heritability
# and add env variance
env1=0.3
env2=0.2

# do both RC and FES, instructions together here
for rep in {1..50}; do
    # trait simulations were performed locally
    time Rscript real-04-simtrait.R --bfile $name -r $rep --sim --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
    time Rscript real-04-simtrait.R --bfile $name -r $rep --sim --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 --fes

    # PCA and LMM runs were actually performed on DCC (see batch.R)
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 --fes
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink --clean

    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 --fes
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean
done

# RC
time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-09-figs.R --bfile $name --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 --final
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2

# FES
time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2
time Rscript real-09-figs.R --bfile $name --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 --final --fes
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2

