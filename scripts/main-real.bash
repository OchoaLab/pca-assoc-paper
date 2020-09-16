# make links to data
# based on: gas-rgls/scripts/data_links.bash
# use LD pruned only

# shared by R script calls further below
name='HoPacAll_ld_prune_1000kb_0.3'
DATA_DIR='/home/viiia/dbs/humanOrigins'

# version for HGDP WGS
name='hgdp_wgs_autosomes_ld_prune_1000kb_0.3'
DATA_DIR='/home/viiia/dbs/hgdp_wgs'

# version for TGP
name='all_phase3_filt-minimal_ld_prune_1000kb_0.3'
DATA_DIR='/home/viiia/dbs/tgp/plink2'

# shared steps (for each in turn)
cd ../data/
mkdir $name
cd $name

ln -s "$DATA_DIR/$name.bed" data.bed
ln -s "$DATA_DIR/$name.bim" data.bim
ln -s "$DATA_DIR/$name.fam" data.fam

# return to scripts dir
cd ../../scripts/
# preprocess with GCTA (makes GRM and max PCs)
time Rscript real-00-preprocess-gcta.R --bfile $name
# 1m25.639s viiiaX6 HO
# 1m26.897s ideapad HGDP
# 19m50.614s ideapad TGP

# get PCs in R, using my formula
# hope for similar average performance, although these PCs are very different than GCTA (except for top 2)
# NOTE: uses ROM version (known bias, does not match popular versions; I think this is best)
time Rscript real-01-pca-test.R --bfile $name
# 6m57.502s viiiaX6 HO
# 5m21.411s ideapad HGDP
# 30m45.184s labbyDuke TGP # ideapad and viiiaX6 ran out of mem (~16G RAM)

# get PCs using plink2
time Rscript real-01-pcs-plink.R --bfile $name
# 0m39.102s ideapad HO
# 0m25.582s ideapad HGDP
# 2m48.309s ideapad TGP

# this creates auxiliary GCTA PCA files (redundant, will delete when done with this analysis)
time Rscript real-02-subset-eigenvec.R --bfile $name
# 0m3.437s ideapad
# same but with standard PCA estimates (from my R code, kinship_std ROM version)
#time Rscript real-02-subset-eigenvec.R --bfile $name --std
# lastly, same but with PCs from plink2
time Rscript real-02-subset-eigenvec.R --bfile $name --plink

# calculates kinship matrix with popkin, to get mean kinship to pass to simtrait
time Rscript real-03-popkin.R --bfile $name
# 4m54.676s ideapad HO
# 8m37.050s ideapad HGDP
# 255m35.418s ideapad TGP

# draws random traits
# 0m2.166s ideapad HO (each rep)
# 0m21.066s ideapad TGP (each rep)
for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep
done

# GCTA runs
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs
    done
done

# GCTA runs
# TGP labbyduke versions
# run reps backwards
for rep in {50..20}; do
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs
    done
done

# # PCA runs (with plink)
# for pcs in {0..90}; do
#     for rep in {1..50}; do
# 	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs
#     done
# done

# PCA runs (with pure plink)
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs --plink
    done
done

# summarizes p-values into AUC and RMSD for each method/rep/pc
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90
# 38+8+5+2+1+3+5+2+1+1+1+1+1+1m ideapad (ran in parts, progressively)
# 61m47.537s + 7m48.435s + 2 + 4 + 10 + 4 (HO PCA complete, GCTA not yet; all 3 machines)
# 654m41.628s + 76m10.185s + 147m53.461s ideapad HGDP
# 33m45.415s ideapad HO + plink pure
# 49m21.674s + 277m34.362s + 56m51.033s ideapad HGDP + plink pure
# ? + 356m44.225s (done) viiiaX6 TGP? PCA pure plink (first parallel run; total user was ? + 1765m9.124s)
# + 146m5.550s ideapad TGP GCTA reps 1-6 (single-threaded due to mem limits)
# + 73m30.572s ideapad TGP GCTA reps 7-9 (single-threaded due to mem limits)
# + 29m11.051s viiiaX6 TGP GCTA reps 10-12 (6 threads)

# read all individual summary tables (tiny files), gather into a master table
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
# 0m25.560s ideapad HO (partial)
# 0m26.445s viiiaX6 HO (partial)
# 0m21.320s labbyduke HO final
# 0m40.056s ideapad HGDP final
# 1m8.339s ideapad HO + plink pure final

# creates final plot for paper!
time Rscript real-09-figs.R --bfile $name
# 0m1.809s ideapad
# for partial runs, to plot only complete reps (best for slowest test: TGP)
time Rscript real-09-figs.R --bfile $name --complete
# compares PCAs only (for internal purposes only)
time Rscript real-09-figs.R --bfile $name --pca

# tests that p-value vectors have the right lengths of m_loci
# to make sure nothing was corrupted due to scripts stopping unexpectedly or incomplete file transfers
# (now test is peformed within real-07-auc-rmsd.R above too, but this can retest everything without recalculating expensive summary statistics).
# --final requires that all files exist!
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90
# 3m12.595s ideapad HO (partial)
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90 --final
# 6m21.249s labbyDuke HO
# 48m44.474s ideapad HGDP
# 13m24.071s ideapad HO + plink pure
# 77m6.361s ideapad HGDP + plink pure

###

# removes redundant, auxiliary GCTA PCA files
time Rscript real-02-subset-eigenvec.R --bfile $name --clean
# 0m0.473s ideapad
#time Rscript real-02-subset-eigenvec.R --bfile $name --clean --std
time Rscript real-02-subset-eigenvec.R --bfile $name --clean --plink

# archive p-values and individual summary files (move out of space that gets synced between computers; these files take up loads of space and are no longer needed)
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90
