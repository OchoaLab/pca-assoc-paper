
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
# link annotations here too
ln -s ../$name/pops-annot.txt pops-annot.txt
# reapply MAF filter, important for consistency of simulated traits, treatment of fixed loci, etc!
time plink2 --bfile data --maf 0.01 --make-bed --out data2
rm data2.log # cleanup
rename data2 data data2.* # replace files


# return to scripts dir
cd ../../scripts/

# preprocess with GCTA (makes GRM and max PCs)
time Rscript real-00-preprocess-gcta.R --bfile $name
# 0m36.076s ideapad HO
# 0m11.685s/1m15.669s viiiaR5 HGDP
# 0m49.143s/5m39.546s viiiaR5 TGP

# get PCs using plink2
time Rscript real-01-pcs-plink.R --bfile $name
# 0m19.751s ideapad HO
# 0m04.303s/0m22.975s viiiaR5 HGDP
# 0m31.934s/3m10.264s viiiaR5 TGP

# calculates kinship matrix with popkin, to get mean kinship to pass to simtrait
time Rscript real-03-popkin.R --bfile $name
# 4m28.038s ideapad HO
# 1m13.881s viiiaR5 HGDP
# 13m20.287s viiiaR5 TGP

# draws random traits
# 0m2.166s ideapad HO (each rep)
# 0m7.475s viiiaR5 HGDP (each rep)
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
# 38m5.081s/314m0.706s HGDP viiiaR5
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
# 39m7.919s/313m59.966s viiiaR5 HGDP
# 46m48.381s/468m21.550s viiiaR5 TGP
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --fes
time Rscript real-09-figs.R --bfile $name --fes
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90 --final --fes
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --fes -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --fes



### TODO NOTE: low-herit and env runs not documented!
# what follows is only gcta-labs model runs

# write labels into a categorical covariates file, to be used by a variant of competitors
time Rscript real-18-make-labs.R --bfile $name

# run GCTA with labels only (no PCs)
for rep in {1..50}; do
    time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs 0 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l
    time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs 0 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l --fes
done

# RC
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-09-figs.R --bfile $name --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 --final -l
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l

# FES
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2
time Rscript real-09-figs.R --bfile $name --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2 -l
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 --final --fes -l
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2 -l -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --herit $h --m_causal_fac $mcf --fes --env1 $env1 --env2 $env2 -l


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

# draw genotypes on DCC, necessary for WGS datasets due to high mem usage
# this is example using "arrays" for replicates
# NOTE: dataset/partition/account/etc are hardcoded here
# sbatch -a 1-50 real-sim.q
# follow up with batch.R runs

# below loop equivalent to real-sim.q/batch.R but run locally
for rep in {1..50}; do
    time Rscript fit-04-draw-geno.R --bfile $name -r $rep --maf_real
    time Rscript real-04-simtrait.R --bfile $name -r $rep --sim
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
# 34m53.484s/335m40.180s viiiaR5 HGDP-sim
# 9m23.831s/88m18.253s viiiaR5 HO-sim
# 46m13.416s/467m16.427s viiiaR5 TGP-sim
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
time Rscript real-09-figs.R --bfile $name
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --final
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90

### FES ###
# below loop equivalent to real-sim.q/batch.R but run locally
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



###################
### KING-cutoff ###
###################

# limited test does only one PC value per method, FES only, real only (no real-sim/tree)

# some constants
pcs_pca=20
pcs_lmm=0
pcs_lmm10=10

#sim=true # change to true/false as needed!
sim=false # change to true/false as needed!

# this handles the biggest difference between real and sim in my scripts...
simflag=''
if [[ "$sim" == true ]]; then simflag='--sim'; fi

# this was saved earlier, ignores "_sim" suffix of last $name
name=$name_king
# some minimal work to get AUCs and RMSDs
# will resimulate traits, meh, though we could have subsetted the original data in theory
time Rscript real-18-make-labs.R --bfile $name
time Rscript real-00-preprocess-gcta.R --bfile $name
time Rscript real-01-pcs-plink.R --bfile $name
time Rscript real-03-popkin.R --bfile $name
for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep --fes
    time Rscript real-04-simtrait.R --bfile $name -r $rep --fes --herit $h --m_causal_fac $mcf
    time Rscript real-04-simtrait.R --bfile $name -r $rep --fes --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
done

# actual runs used batch.R, but it's equivalent to the below steps
# a bit overkill making all PCs, but faster than rewriting the code for this special case
time Rscript real-02-subset-eigenvec.R --bfile $name --plink
for rep in {1..50}; do
    time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs_pca --fes
    time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs_pca --fes --herit $h --m_causal_fac $mcf
    time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs_pca --fes --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
done
time Rscript real-02-subset-eigenvec.R --bfile $name --clean --plink
time Rscript real-02-subset-eigenvec.R --bfile $name
for rep in {1..50}; do
    time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs_lmm --fes
    time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs_lmm --fes --herit $h --m_causal_fac $mcf
    time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs_lmm --fes --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
    time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs_lmm10 --fes --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
    time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs_lmm --fes --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l
done
time Rscript real-02-subset-eigenvec.R --bfile $name --clean


# FES h=0.8
# gather results into a single table (again overkill but meh)
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs $pcs_pca --fes
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs $pcs_pca --fes
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs $pcs_pca --fes # NOTE not --final
# note extra flags for partial data
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_pca --fes -s -m pca
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_lmm --fes -s -m lmm

# FES h=0.3
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs $pcs_pca --fes --herit $h --m_causal_fac $mcf
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs $pcs_pca --fes --herit $h --m_causal_fac $mcf
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs $pcs_pca --fes  --herit $h --m_causal_fac $mcf # NOTE not --final
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_pca --fes -s -m pca --herit $h --m_causal_fac $mcf
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_lmm --fes -s -m lmm --herit $h --m_causal_fac $mcf

# FES h=0.3 env
# in these cases there is an assumption that $pcs_pca is the highest of all values, which is true but beware!
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs $pcs_pca --fes --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs $pcs_pca --fes --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs $pcs_pca --fes  --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l # NOTE not --final
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_pca --fes -s -m pca --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_lmm --fes -s -m lmm --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_lmm10 --fes -s -m lmm --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs $pcs_lmm --fes -s -m lmm-labs --herit $h --m_causal_fac $mcf --env1 $env1 --env2 $env2 -l
