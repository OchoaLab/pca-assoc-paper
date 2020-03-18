library(optparse)
library(genio)       # to write BED files for external software
library(popkinsuppl) # for PCA's kinship estimator

#setwd("../scripts")
setwd("/dscrhome/yy222/gas-rgls-master/scripts1")

# standard code for a complex trait and an admixed population
source('sim_geno_trait_k3.R')

# load new functions from external scripts
source('kinship_to_evd.R')
source('gas_lm_optim.R')
source('gas_pca_optim.R')
source('gas_lmm_gcta.R')
source('paths.R')
source("gas_plots.R")

# number of replicates per run
rep <- 2
# number of PCs to explore
n_pcs_max <- 90
# alternate path for GCTA binary
gcta_bin <- "/dscrhome/yy222/gcta_1.92.4beta1/gcta64"

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-n", "--n_ind"), type = "integer", default = 1000, # CHANGE: or 100
                help = "number of individuals", metavar = "int"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 100000, 
                help = "number of loci", metavar = "int"),
    make_option(c("-k", "--k_subpops"), type = "integer", default = 10, 
                help = "admixture intermediate subpopulations", metavar = "int"),
    make_option(c("-f", "--fst"), type = "double", default = 0.1, 
                help = "FST (fixation index)", metavar = "double"),
    make_option(c("--bias_coeff"), type = "double", default = 0.5, 
                help = "admixture bias coeff", metavar = "double"),
    make_option(c("-g", "--generations"), type = "integer", default = 1, 
                help = "number of generations, for realistic local kinship", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal", type = "integer", default = 100, # CHANGE: or 10
                help = "num causal loci", metavar = "int"),
    make_option(c("-t", "--threads"), type = "integer", default = 1, 
                help = "number of threads (affects GCTA only)", metavar = "int"),
    make_option("--debug", action = "store_true", default = FALSE, 
                help = "debug mode (GCTA is fully verbose)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_ind <- opt$n_ind
m_loci <- opt$m_loci
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
generations <- opt$generations
m_causal <- opt$m_causal
herit <- opt$herit
threads <- opt$threads
debug <- opt$debug

# move now to where we want outputs to be
## dir.create('test_yiqi')
## setwd('test_yiqi')
setwd( paste0("/dscrhome/yy222/LMM_PCA_rep_n_", n_ind, "/rep1") )

# output path for BED files and figure
name_out <- paste0(
    'gas',
    '-n', n_ind,
    '-m', m_loci,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-mc', m_causal,
    '-h', herit,
    '-t', threads,
    '-g', generations
)

# in all cases go up to n_pcs_max + 1 to include the zero PC cases
M_rmsd <- matrix(0, rep, n_pcs_max + 1)
M_auc <- matrix(0, rep, n_pcs_max + 1)
M_rmsd_gcta <- matrix(0, rep, n_pcs_max + 1)
M_auc_gcta <- matrix(0, rep, n_pcs_max + 1)

for (i in 1 : rep){
    message('rep: ', i)

    ############
    ### SIMS ###
    ############
    
    # simulate trait!
    # can be done inside an inner loop to simulate thousands of traits per genotype dataset
    message('simtrait')
    # NOTE: using kinship version, only choice for real data
    obj <- sim_geno_trait_k3(
        n_ind = n_ind,
        m_loci = m_loci,
        m_causal = m_causal,
        k_subpops = k_subpops,
        bias_coeff = bias_coeff,
        generations = generations,
        herit = herit,
        verbose = TRUE,
        fst = fst
    )
    X <- obj$X
    trait <- obj$trait
    causal_indexes <- obj$causal_indexes
    
    # write plink data for GCTA
    plink_data <- write_plink(name_out, X, pheno = trait)
    write_phen(name_out, plink_data$fam)

    # estimate kinship (std) and eigenvectors for PCA
    kinship_estimate_old <- kinship_std(X) # estimate kinship the old way
    eigenvectors_estimate_old <- kinship_to_evd( kinship_estimate_old ) # get all eigenvalues

    ########################
    ### ASSOCIATION TEST ###
    ########################

    message("LM")
    
    obj <- gas_lm_optim(X, trait)
    # NOTE zero PCS is position 1
    M_rmsd[i, 1] <- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
    M_auc[i, 1] <- pvals_to_pr_auc(obj$pvals, causal_indexes)

    # GCTA with zero PCs
    message("GCTA")
    # make GRM (shared by all GCTA runs)
    obj <- gas_lmm_gcta_kin(gcta_bin, name_out, debug = debug)
    # actual GWAS run
    obj <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads, debug = debug)
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
        gas_lmm_gcta_pca(gcta_bin, name_out, n_pcs = pcs, debug = debug)

        # actual GWAS
        obj <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads, debug = debug, n_pcs = pcs)
        # NOTE r PCS is position r+1
        M_rmsd_gcta[i, pcs + 1] <- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
        M_auc_gcta[i, pcs + 1] <- pvals_to_pr_auc(obj$pvals, causal_indexes)
        
        # clean up when we're done with gcta
        delete_files_gcta(name_out) # deletes GAS table only, not GRM (as desired)
        # these get deleted separately
        delete_files_gcta_pca(name_out, n_pcs = pcs)
    }

    # final cleanup
    # delete more files after entire GCTA analysis is complete
    # delete GRM 
    delete_files_gcta(name_out, grm = TRUE)
    # delete plink files
    delete_files_plink(name_out)
    # delete phen file too
    file_phen <- paste0(name_out, '.phen')
    # test that it's actually there!
    if( file.exists(file_phen) ) {
        # remove if it was there
        invisible( file.remove(file_phen) )
    } else {
        warning('File to remove did not exist: ', file_phen)
    }
}

# write outputs
write.table(M_auc_gcta, file = paste0("auc_gcta_n_", n_ind, ".txt") )
write.table(M_rmsd_gcta, file = paste0("rmsd_gcta_n_", n_ind, ".txt") )
write.table(M_auc, file = paste0("auc_pca_n_", n_ind, ".txt") )
write.table(M_rmsd, file = paste0("rmsd_pca_n_", n_ind, ".txt") )
