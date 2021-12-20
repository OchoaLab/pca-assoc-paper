
# creates a table other scripts might use, especially when creating figures, for consistent human-facing naming, and loops across all datasets
# writes data/datasets.txt
Rscript real-00-datasets.R

# make links to data
# based on: gas-rgls/scripts/data_links.bash
# use LD pruned only

# shared by R script calls further below
name='HoPacAll_ld_prune_1000kb_0.3_maf-0.01'
DATA_DIR='/home/viiia/dbs/humanOrigins'

# version for HGDP WGS
name='hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01'
DATA_DIR='/home/viiia/dbs/hgdp_wgs'

# version for TGP
name='tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01'
DATA_DIR='/home/viiia/dbs/tgp-nygc'

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

# NOTE: these actually run on DCC through batch.R (requires manual edits)

# redo if needed
time Rscript real-02-subset-eigenvec.R --bfile $name --dcc --plink
# PCA runs (with pure plink)
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs --fes
    done
done
time Rscript real-02-subset-eigenvec.R --bfile $name --dcc --clean --plink

# redo if needed
time Rscript real-02-subset-eigenvec.R --bfile $name --dcc
# GCTA runs
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs --fes
    done
done
time Rscript real-02-subset-eigenvec.R --bfile $name --dcc --clean

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

# what follows has been run for TGP only!

# manually copied annotations file from Storey Lab project
# for TGP
cp ~/docs/ochoalab/fst-human/data/tgp-nygc-autosomes/pops-annot.txt ../data/$name/
# for HGDP
cp ~/docs/ochoalab/fst-human/data/hgdp_wgs_autosomes/pops-annot.txt ../data/$name/
# for HO
cp ~/docs/ochoalab/fst-human/data/HoPacAll/pops-annot.txt ../data/$name/

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
sbatch -a 1-50 real-sim.q

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


###############
### GLOBALS ###
###############

# cross-dataset analyses, including all of: sim, real, and real-sim
# this sequence also creates all of the main figures/tables in order of appearance

# reports actual dimensions (n_ind, m_loci, K, and m_causal) used in all datasets
# (validates every replicate too! for m_causal compares both FES and RC!)
# writes data/dimensions.txt
time Rscript real-14-dimensions.R
# 0m24.674s viiiaR5

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
# Inflation factor: pca-plink-pure: 2.96003533785624
# Inflation factor: gcta: 0.87992175683473
# Inflation factor: gcta: 0.984541944592377

# main statistical evaluations between methods
time Rscript real-13-stats.R

# final plot gathers all three datasets into a single multipanel figure
time Rscript real-15-plots-big.R
time Rscript real-15-plots-big.R --fes
time Rscript real-15-plots-big.R --real
time Rscript real-15-plots-big.R --real --fes
time Rscript real-15-plots-big.R --real_sim
time Rscript real-15-plots-big.R --real_sim --fes

# calculages popkin and eigensoft eigenvectors, calculates TW stats, makes plot
time Rscript all-02-eigen.R
# 30m14.103s viiiaR5 first time

# a comparison of RMSD and lambda across ALL datasets (including FES and RC traits)
# this version that fits top half only (makes most sense for our goal of talking mostly about inflation)
time Rscript real-11-inflation-across-datasets.R
# model fit:
# rmsd ~ a * (lambda^b - 1) / (lambda^b + 1)
#         a         b 
# 0.5629478 0.6208654 
# threshold map (sigmoidal): lambda = 1.05, RMSD = 0.00853
# Inverse threshold map (sigmoidal): RMSD = 0.01, lambda = 1.06
# log-linear approx: log(lambda) = RMSD * 5.72
# threshold map (log-linear): lambda = 1.05, RMSD = 0.00853
