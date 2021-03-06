
# shared steps (for each real dataset in turn)
cd ../data/
mkdir $name
cd $name

ln -s "$DATA_DIR/$name.bed" data.bed
ln -s "$DATA_DIR/$name.bim" data.bim
ln -s "$DATA_DIR/$name.fam" data.fam
ln -s "$DATA_DIR/pops-annot.txt" pops-annot.txt

# estimate local kinship with KING-robust
time plink2 --bfile data --make-king triangle bin4 --out data
# 0m2.411s/0m16.069s HGDP
# 0m3.979s/0m32.770s HO
# 0m15.863s/2m4.843s TGP
# cleanup
rm data.log

# now create kinship-filtered data
time plink2 --bfile data --king-cutoff data 0.02209709 --make-bed --out data-king-cutoff
# 0m2.242s HGDP
# 0m1.306s HO
# 0m4.910s TGP
# cleanup
rm data-king-cutoff.king.cutoff.{in,out}.id data-king-cutoff.log
# make new folder with this filtered data
name_king=$name'_king-cutoff-4'
mkdir ../$name_king
mv data-king-cutoff.bed ../$name_king/data.bed
mv data-king-cutoff.bim ../$name_king/data.bim
mv data-king-cutoff.fam ../$name_king/data.fam
cd ../$name_king
# reapply MAF filter, important for consistency of simulated traits, treatment of fixed loci, etc!
time plink2 --bfile data --maf 0.01 --make-bed --out data2
rm data2.log # cleanup
rename data2 data data2.* # replace files


# return to scripts dir
cd ../../scripts/

# preprocess with GCTA (makes GRM and max PCs)
time Rscript real-00-preprocess-gcta.R --bfile $name
# 0m36.076s ideapad HO
# 6m17.134s/66m7.781s viiiaR5 HGDP
# 0m49.143s/5m39.546s viiiaR5 TGP

# get PCs using plink2
time Rscript real-01-pcs-plink.R --bfile $name
# 0m19.751s ideapad HO
# 0m9.476s/0m47.619s viiiaR5 HGDP
# 0m31.934s/3m10.264s viiiaR5 TGP

# calculates kinship matrix with popkin, to get mean kinship to pass to simtrait
time Rscript real-03-popkin.R --bfile $name
# 4m28.038s ideapad HO
# 2m25.247s viiiaR5 HGDP
# 13m20.287s viiiaR5 TGP

# draws random traits
# 0m2.166s ideapad HO (each rep)
# 0m9.504s viiiaR5 HGDP (each rep)
# 0m18.936s viiiaR5 TGP (each rep)
for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep
done

# gather minimal data to send to cluster, including FES versions (not yet made in this script)
# (tar -h: "Follow symlinks; archive and dump the files they point to.")
# cd ../data/$name/
# tar -chzf data.tgz data.{bed,bim,fam} data.grm.* data-n_pcs_90.eigenvec data-plink-n_pcs_90.eigenvec rep-*/data.phen rep-*/fes/data.phen 

# create auxiliary PCA files from plink2 (redundant, will delete when done with this analysis)
time Rscript real-02-subset-eigenvec.R --bfile $name --plink
# PCA runs (with pure plink)
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs
    done
done
# removes redundant, auxiliary plink2 PCA files
time Rscript real-02-subset-eigenvec.R --bfile $name --clean --plink

# create auxiliary GCTA PCA files (redundant, will delete when done with this analysis)
time Rscript real-02-subset-eigenvec.R --bfile $name
# 0m3.437s ideapad
# GCTA runs
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs
    done
done
# removes redundant, auxiliary GCTA PCA files
time Rscript real-02-subset-eigenvec.R --bfile $name --clean
# 0m0.473s ideapad

# summarizes p-values into AUC and RMSD for each method/rep/pc
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90
# all are plink + GCTA
# 25m13.604s HO viiiaX6 6 threads
# 39m51.574s/407m5.064s HGDP viiiaR5
# 48m4.522s/460m21.524s TGP viiiaR5

# read all individual summary tables (tiny files), gather into a master table
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
# 0m42.030s viiiaX6 HO
# 0m49.292s viiiaR5 HGDP
# 0m51.765s viiiaR5 TGP

# creates final plot for paper!
time Rscript real-09-figs.R --bfile $name
# 0m1.809s ideapad
# # OBSOLETE for partial runs, to plot only complete reps (best for slowest test: TGP)
# time Rscript real-09-figs.R --bfile $name --complete

# tests that p-value vectors have the right lengths of m_loci
# to make sure nothing was corrupted due to scripts stopping unexpectedly or incomplete file transfers
# (now test is peformed within real-07-auc-rmsd.R above too, but this can retest everything without recalculating expensive summary statistics).
# --final requires that all files exist!
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90 --final
# 7m57.939s viiiaX6 HO
# 9m17.899s viiiaR5 HGDP
# 12m20.459s viiiaR5 TGP

# archive p-values and individual summary files (move out of space that gets synced between computers; these files take up loads of space and are no longer needed)
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90

# calculate and store MAF distributions, useful not just for a plot but also for simulating data from real datasets
time Rscript real-16-mafs.R --bfile $name

###########
### FES ###
###########

# addition to complement existing analysis with alternative trait from FES model
# not sure if it will make a meaningful difference
# NOTE: many steps that depend on genotypes only aren't redone (are shared from prev run)

for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep --fes
done

# NOTE: these were actually run on cluster through batch.R (requires manual edits)

# redo if needed
time Rscript real-02-subset-eigenvec.R --bfile $name --plink
# PCA runs (with pure plink)
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs --fes
    done
done
time Rscript real-02-subset-eigenvec.R --bfile $name --clean --plink

# redo if needed
time Rscript real-02-subset-eigenvec.R --bfile $name
# GCTA runs
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs --fes
    done
done
time Rscript real-02-subset-eigenvec.R --bfile $name --clean

time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90 --fes
# 56m22.873s/342m28.028s viiiaR5 HGDP
# 46m48.381s/468m21.550s viiiaR5 TGP
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --fes
time Rscript real-09-figs.R --bfile $name --fes
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90 --final --fes
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --fes -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --fes


###################
### FIT FOR SIM ###
###################

# we already had popkin estimate because we needed params for simtrait (all real datasets)

# this script calculates subpopulations kinship matrix, which we'll fit tree to
# also produces a nice plot that validates the estimate
time Rscript fit-01-popkin-subpops.R --bfile $name
# 4s ideapad
# NOTE: label/margin sizes tuned for TGP, terrible for HGDP and HO

# fit tree to real data, visualize success
time Rscript fit-02-tree.R --bfile $name

# set up simulation params Q and Psi tree
# also calculates FST for paper table and kinship_mean for undiff_af
time Rscript fit-03-sim-pop.R --bfile $name

# for extra convenience, make a "copy" of the real MAF file into the simulations dir (better for cluster work)
ln -s ../$name/maf.RData ../data/$name"_sim"/maf-real.RData 

# from here on, the real data name needs a '_sim' suffix
# the rest of the pipeline is the same as for the regular simulations
name=$name"_sim"

# # OPTIONAL experiments to help me decide how to fit this all
# # can skip: ultimately does not produce required outputs, only influenced choices hardcoded into fit-04-draw-geno.R now
# time Rscript fit-06-bigger-fitting.R --bfile $name
# # 2m23.877s viiiaR5 HGDP
# # 6m25.336s viiiaR5 HO
# time Rscript fit-06-bigger-fitting.R --bfile $name --beta
# # 3m23.708s viiiaR5 HGDP
# # 20m50.410s viiiaR5 HO
# time Rscript fit-06-bigger-fitting.R --bfile $name --distr
# # 6m54.920s viiiaR5 TGP
# # 1m41.267s viiiaR5 HGDP
# # 6m45.279s viiiaR5 HO
# time Rscript fit-06-bigger-fitting.R --bfile $name --kin
# # 2m3.379s viiiaR5 HGDP
# # 8m2.776s viiiaR5 TGP
# # 8m9.043s viiiaR5 HO
# # and plotter code (after all three are done)
# time Rscript fit-07-bigger-fitting-plot.R --bfile $name
# # runs with final m (default is smaller) to have more accurate error measurements (way slower though)
# # didn't run for HO (all previous signs suggested it wasn't needed, though ultimately these weren't needed either)
# time Rscript fit-06-bigger-fitting.R --bfile $name --kin -m 100000
# # 19m41.688s viiiaR5 HGDP
# # 72m39.955s viiiaR5 TGP
# time Rscript fit-07-bigger-fitting-plot.R --bfile $name -m 100000

# draw genotypes on DCC
# this is example using "arrays" for replicates
# NOTE: dataset/partition/account/etc are hardcoded here
# sbatch -a 1-50 real-sim.q

# actually simulate data
for rep in {1..50}; do
    time Rscript fit-04-draw-geno.R --bfile $name -r $rep --maf_real
    time Rscript sim-02-sim-trait.R --bfile $name -r $rep
    time Rscript real-00-preprocess-gcta.R --bfile $name/rep-$rep
    time Rscript real-01-pcs-plink.R --bfile $name/rep-$rep
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink --clean
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean
done
# new steps to make sure simulation is as expected (validates not just tree but also Q matrix and their respective alignment)
# estimates from rep-1 only!
time Rscript real-03-popkin.R --bfile $name --sim
# an FST validation (model vs sim)
time Rscript fit-05-fst-validate.R --bfile $name
# make plot based on those estimates
time Rscript fit-10-plot-real-vs-sim.R --bfile $name
# get MAFs for comparison plot too (only needs rep-1)
time Rscript real-16-mafs.R --bfile $name --sim

# estimate local kinship with KING-robust, rep-1 only
cd ../data/$name/rep-1
time plink2 --bfile data --make-king triangle bin4 --out data
# cleanup
rm data.log 
# return to scripts dir
cd ../../../scripts/

time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90
# 7m28.614s/38m52.149s viiiaR5 12 threads HGDP-sim
# 39m0.437s/408m46.361s viiiaR5 HGDP-sim full-m
# 9m23.831s/88m18.253s viiiaR5 HO-sim full-m
# 46m13.416s/467m16.427s viiiaR5 TGP-sim full-m
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
time Rscript real-09-figs.R --bfile $name
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --final
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90

### FES ###
for rep in {1..50}; do
    time Rscript sim-02-sim-trait.R --bfile $name -r $rep --fes
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



###################
### KING-cutoff ###
###################

# limited test does only one PC value per method, FES only, real only (no real-sim/tree)

# some constants
pcs_pca=20
pcs_lmm=0

# this was saved earlier, ignores "_sim" suffix of last $name
name=$name_king
# some minimal work to get AUCs and RMSDs
# will resimulate traits, meh, though we could have subsetted the original data in theory
time Rscript real-00-preprocess-gcta.R --bfile $name
time Rscript real-01-pcs-plink.R --bfile $name
time Rscript real-03-popkin.R --bfile $name
for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep --fes
done
# a bit overkill making all PCs, but faster than rewriting the code for this special case
time Rscript real-02-subset-eigenvec.R --bfile $name --plink
for rep in {1..50}; do
    time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs_pca --fes
done
time Rscript real-02-subset-eigenvec.R --bfile $name --clean --plink
# no PCS for GCTA here only
for rep in {1..50}; do
    time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs_lmm --fes
done
# gather results into a single table (again overkill but meh)
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs $pcs_pca --fes
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs $pcs_pca --fes
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs $pcs_pca --fes # NOTE not --final
# note extra flags for partial data
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_pca --fes -s -m pca -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_pca --fes -s -m pca
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_lmm --fes -s -m lmm -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_lmm --fes -s -m lmm

