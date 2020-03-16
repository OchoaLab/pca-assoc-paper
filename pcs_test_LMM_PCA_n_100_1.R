library(optparse)

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-p", "--plot"), action = "store_true", default = FALSE, 
                help = "create plot of existing data"),
    make_option("--simple", action = "store_true", default = FALSE, 
                help = "simple plot version (for a grant proposal)"),
    make_option(c("-l", "--lite"), action = "store_true", default = FALSE, 
                help = "run lite version (focuses on best and least redundant methods)"),
    make_option(c("-n", "--n_ind"), type = "integer", default = 100, 
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
    make_option("--m_causal", type = "integer", default = 10, 
                help = "num causal loci", metavar = "int"),
    make_option(c("-t", "--threads"), type = "integer", default = 1, 
                help = "number of threads (affects GCTA only)", metavar = "int"),
    make_option("--debug", action = "store_true", default = FALSE, 
                help = "debug mode (GCTA and GEMMA are fully verbose)")
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
lite <- opt$lite
threads <- opt$threads
simple <- opt$simple
debug <- opt$debug

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


 library(readr)       # to write kinship matrix
    library(tibble)      # to store data
    library(genio)       # to write BED files for external software
    library(popkin)      # to estimate kinship in RGLS
    library(popkinsuppl) # for PCA's kinship estimator
    library(lfa)         # GWAS gcatest
    library(gcatest)     # GWAS gcatest
    #library(qvalue)      # multiple hypothesis tests

setwd("/dscrhome/yy222/gas-rgls-master/scripts1")

    # standard code for a complex trait and an admixed population
    source('sim_geno_trait_k3.R')

    # load new functions from external scripts
    source('kinship_to_evd.R')
    source('gas_lm_optim.R')
    source('gas_pca_optim.R')
    source('gas_rgls.R')
    source('gas_lmm_gemma.R')
    source('gas_lmm_emmax.R')
    source('gas_lmm_gcta.R')
    source('paths.R')
        source("gas_plots.R")
# number of replicates per run
rep <- 2
# number of PCs to explore
n_pcs_max <- 2
# alternate path for GCTA binary
gcta_bin<-"/dscrhome/yy222/gcta_1.92.4beta1/gcta64"


############
### SIMS ###
############
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
    admix_proportions <- obj$admix_proportions
    kinship <- obj$kinship
    trait <- obj$trait
    causal_indexes <- obj$causal_indexes
    
    ## how we'd get these values if we weren't setting them on command line
    ## n_ind <- ncol(X)
    ## m_loci <- nrow(X)
    ## k_subpops <- ncol(admix_proportions)
    ## sanity checks instead
    stopifnot( ncol(X) == n_ind )
    stopifnot( nrow(X) == m_loci )
    stopifnot( ncol(admix_proportions) == k_subpops )
    
    # write BED version for external code
    plink_data <- write_plink(name_out, X, pheno = trait, verbose = FALSE)

    # write phenotype file
    # shared by EMMAX and GCTA
    write_phen(name_out, plink_data$fam, verbose = FALSE)# in all cases go up to n_pcs_max + 1 to include the zero PC cases

M_rmsd <- matrix(0, rep, n_pcs_max + 1)
M_auc <- matrix(0, rep, n_pcs_max + 1)
M_rmsd_gcta <- matrix(0, rep, n_pcs_max + 1)
M_auc_gcta <- matrix(0, rep, n_pcs_max + 1)
eigenvectors <- kinship_to_evd(kinship)

for (i in 1 : rep){
    message('rep: ', i)
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
    admix_proportions <- obj$admix_proportions
    kinship <- obj$kinship
    trait <- obj$trait
    causal_indexes <- obj$causal_indexes    # data dimensions
    # simulation is only BEDMatrix, but meh
    m_loci <- if (class(X) == 'BEDMatrix') ncol(X) else nrow(X)
    

     ## how we'd get these values if we weren't setting them on command line
    ## n_ind <- ncol(X)
    ## m_loci <- nrow(X)
    ## k_subpops <- ncol(admix_proportions)
    ## sanity checks instead
    stopifnot( ncol(X) == n_ind )
    stopifnot( nrow(X) == m_loci )
    stopifnot( ncol(admix_proportions) == k_subpops )

   plink_data <- write_plink(name_out, X, pheno = trait, verbose = FALSE)

  write_phen(name_out, plink_data$fam, verbose = FALSE)

        eigenvectors <- kinship_to_evd(kinship)

        # write true kinship matrix in correct format
        # shared by at least GEMMA and EMMAX
        file_true_kin <- paste0(name_out, '.true_kinship.txt')
        kinship_tibble <- as_tibble( 2 * kinship, .name_repair = 'minimal')
        write_tsv(kinship_tibble, file_true_kin, col_names = FALSE)
    
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
    obj <- gas_lmm_gcta_kin(gcta_bin, name_out, debug = debug)
    obj <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads, debug = debug)
    # NOTE zero PCS is position 1
    M_rmsd_gcta[i, 1] <- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
    M_auc_gcta[i, 1] <- pvals_to_pr_auc(obj$pvals, causal_indexes)
    
    for (pcs in 1 : n_pcs_max){
        message('PCs: ', pcs)
        
        message("PCA")
        indexes <- 1 : pcs
        kinship_estimate_old <- kinship_std(X) # estimate kinship the old way
        eigenvectors_estimate_old <- kinship_to_evd( kinship_estimate_old ) # get all eigenvalues
        obj <- gas_pca_optim(X, trait, eigenvectors_estimate_old[, indexes ])
	#    file_gcta_kin <- obj$file	
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
}


setwd("/dscrhome/yy222/LMM_PCA_rep_n_100/rep1")
write.table(M_auc_gcta, file = "auc_gcta_n_100.txt")
write.table(M_rmsd_gcta, file = "rmsd_gcta_n_100.txt")
write.table(M_auc, file = "auc_pca_n_100.txt")
write.table(M_rmsd, file = "rmsd_pca_n_100.txt")