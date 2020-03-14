# for terminal options
library(optparse)
library(BEDMatrix) # to load real genotypes with low memory usage
library(simtrait)  # to simulate complex trait
library(genio)     # to write BED files for external software

# load new functions from external scripts
setwd("/dscrhome/yy222/gas-rgls-master/scripts1")
#setwd("../scripts")
source('gas_lm_optim.R')
source('gas_pca_optim.R')
source('gas_lmm_gcta.R')
source('paths.R')
source('gas_plots.R')

# number of replicates per run
rep <- 2
# number of PCs to explore
n_pcs_max <- 2 # 90
# alternate path for GCTA binary
gcta_bin<-"/hpc/group/biostat/yy222/gcta_1.92.4beta1/gcta64"

############
### ARGV ###
############
# define options
option_list = list(
    make_option(c("-i", "--input"), type = 'character', default = "human_origins_and_pacific_public", 
                help = "Path for input genotypes (BED/BIM/FAM files, excluding extension)"),
    make_option("--debug", action = "store_true", default = FALSE, 
                help = "debug mode (GCTA is fully verbose)"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal", type = "integer", default = 100, 
                help = "num causal loci", metavar = "int"),
    make_option(c("-t", "--threads"), type = "integer", default = 1, 
                help = "number of threads (affects GCTA only)", metavar = "int")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name_in <- opt$input
m_causal <- opt$m_causal
herit <- opt$herit
threads <- opt$threads
debug <- opt$debug

############
### SIMS ###
############
setwd("/dscrhome/yy222/human-differentiation-manuscript-master/data")
#setwd('../../human-differentiation-manuscript/data/')

# load precomputed data: eigenvectors_estimate_old, kinship_estimate
file_kinship_data <- paste0(name_in, '.RData')
message('loading: ', file_kinship_data)
load(file_kinship_data)

# load genotypes
message('BEDMatrix')
X <- BEDMatrix(name_in)

# in all cases, fam$fam are good labels
# load data with genio
fam <- read_fam(name_in)

# in all cases go up to n_pcs_max + 1 to include the zero PC cases
M_rmsd <- matrix(0, rep, n_pcs_max + 1)
M_auc <- matrix(0, rep, n_pcs_max + 1)
M_rmsd_gcta <- matrix(0, rep, n_pcs_max + 1)
M_auc_gcta <- matrix(0, rep, n_pcs_max + 1)

for (i in 1 : rep){
    message('rep: ', i)
    # simulate trait!
    # can be done inside an inner loop to simulate thousands of traits per genotype dataset
    message('simtrait')
    # NOTE: using kinship version, only choice for real data
    obj_trait <- sim_trait(
        X = X,
        m_causal = m_causal,
        herit = herit,
        kinship = kinship_estimate
    )
    trait <- obj_trait$trait
    causal_indexes <- obj_trait$causal_indexes
    
    # data dimensions
    # simulation is only BEDMatrix, but meh
    m_loci <- if (class(X) == 'BEDMatrix') ncol(X) else nrow(X)
    
    # add phenotype to FAM table
    # safe since initial phenotypes were all 1 (essentially no data)
    fam$pheno <- trait
    # write phenotype file
    # for GCTA
    write_phen(name_in, fam)
    
    message("LM")
    
    obj <- gas_lm_optim(X, trait)
    # NOTE zero PCS is position 1
    M_rmsd[i, 1] <- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
    M_auc[i, 1] <- pvals_to_pr_auc(obj$pvals, causal_indexes)

    # GCTA with zero PCs
    message("GCTA")
    
    obj <- gas_lmm_gcta(gcta_bin, name_in, m_loci = m_loci, threads = threads, debug = debug)
    # NOTE zero PCS is position 1
    M_rmsd_gcta[i, 1] <- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
    M_auc_gcta[i, 1] <- pvals_to_pr_auc(obj$pvals, causal_indexes)
    
    for (pcs in 1 : n_pcs_max){
        message('PCs: ', pcs)
        
        message("PCA")
        indexes <- 1 : pcs
        obj <- gas_pca_optim(X, trait, eigenvectors_estimate_old[, indexes ])
        # NOTE r PCS is position r+1
        M_rmsd[i, pcs + 1] <- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
        M_auc[i, pcs + 1] <- pvals_to_pr_auc(obj$pvals, causal_indexes)
        
        # GCTA with PCA
        message("GCTAp")
        
        # gets PCs, saves them on a file
        gas_lmm_gcta_pca(gcta_bin, name_in, n_pcs = pcs, debug = debug)

        # actual GWAS
        obj <- gas_lmm_gcta(gcta_bin, name_in, m_loci = m_loci, threads = threads, debug = debug, n_pcs = pcs)
        # NOTE r PCS is position r+1
        M_rmsd_gcta[i, pcs + 1] <- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
        M_auc_gcta[i, pcs + 1] <- pvals_to_pr_auc(obj$pvals, causal_indexes)
        
        # clean up when we're done with gcta
        delete_files_gcta(name_in) # deletes GAS table only, not GRM (as desired)
        # these get deleted separately
        delete_files_gcta_pca(name_in, n_pcs = pcs)
    }
}

write.table(M_auc_gcta, file = "auc_gcta.txt")
write.table(M_rmsd_gcta, file = "rmsd_gcta.txt")
write.table(M_auc, file = "auc_pca.txt")
write.table(M_rmsd, file = "rmsd_pca.txt")
