# for terminal options
library(optparse)
setwd("/dscrhome/yy222/gas-rgls-master/scripts")
# should be a command line option, meh
verbose <- TRUE

############
### ARGV ###
############

# define options
option_list = list(
  make_option(c("-i", "--input"), type = 'character', default = NA, 
              help = "Path for input genotypes (BED/BIM/FAM files, excluding extension)"),
  make_option(c("-p", "--plot"), action = "store_true", default = FALSE, 
              help = "create plot of existing data"),
  make_option("--simple", action = "store_true", default = FALSE, 
              help = "simple plot version (for a grant proposal)"),
  make_option(c("-l", "--lite"), action = "store_true", default = FALSE, 
              help = "run lite version (focuses on best and least redundant methods)"),
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
name_in <- "/dscrhome/yy222/gas-rgls-master/scripts/human_origins_and_pacific_public"

if (is.na(name_in))
  stop('`-i/--input` is mandatory!')
m_causal <- opt$m_causal
herit <- opt$herit
lite <- opt$lite
threads <- opt$threads
simple <- opt$simple

# output path for BED files and figure
# NOTE: here name_in may (and does) have slashes because it's in a directory, so this only works if gas-* is not a prefix as in the original script
name_out <- paste0(
  name_in,
  '-gas',
  '-mc', m_causal,
  '-h', herit,
  '-t', threads
)
# mark simplified plots too
# NOTE: again use suffixes only, name_out has directory slashes
if (lite)
  name_out <- paste0(name_out, '-lite')




library(readr) # to read style



# plotting code
source('gas_plots.R')

# load style table (maps methods to colors and line styles)
style <- read_tsv('style.txt', col_types = 'cci')



# else run simulation

library(BEDMatrix) # to load real genotypes with low memory usage
library(simtrait)  # to simulate complex trait
library(readr)     # to write kinship matrix
library(tibble)    # to store data
library(genio)     # to write BED files for external software
library(gcatest)   # GWAS gcatest

# load new functions from external scripts
source('gas_lm_optim.R')
source('gas_pca_optim.R')
source('gas_rgls.R')
source('gas_lmm_gemma.R')
source('gas_lmm_emmax.R')
source('gas_lmm_gcta.R')
source('paths.R')

############
### SIMS ###
############


load("/dscrhome/yy222/real_data_set/file_kinship_data_human.RData")



# load genotypes
if (verbose)
  message('BEDMatrix')
X <- BEDMatrix(name_in)

# in all cases, fam$fam are good labels
# load data with genio
fam <- read_fam(name_in, verbose = verbose)
M_rmsd<-matrix(0,5,100)
M_auc<-matrix(0,5,100)

style <- read_tsv('style.txt', col_types = 'cci')

  for (j in 1:5){  # simulate genotypes and trait as usual
    
    # simulate trait!
    # can be done inside an inner loop to simulate thousands of traits per genotype dataset
    if (verbose)
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
    times_gas <- tibble(.rows = 1)
    pvals <- tibble(.rows = m_loci)
    
    ########################
    ### ASSOCIATION TEST ###
    ########################
    for (pcs in 2:30){
    
    # since true K isn't known, only K=10 is tested
    
    name <- "PCA"
    message(name)
    
    indexes <- 1 : 10
    obj <- gas_pca_optim(X, trait, eigenvectors_estimate_old[, pcs ])
    times_gas[[name]] <- obj$runtime
    pvals[[name]] <- obj$pvals
    
    
    
    
    style <- read_tsv('style.txt', col_types = 'cci')
    
    
    
    
    
    pvals_auc_collect  <- function (name_out, causal_indexes, pval_tab, times_kin, times_gas, style, cex_leg = 0.7) {
      # we shall generate other metrics in the plots below
      pval_rmsd <- list() # p-value mean distance from desired (under null)
      pr_auc <- list() # AUCs for P-R curves
      
      # reorder input tables as desired (style table)
      # also performs sanity checks (checks if style table is missing names)
      pval_tab <- reorder_tibble_names(style$name, pval_tab)
      times_kin <- reorder_tibble_names(style$name, times_kin)
      times_gas <- reorder_tibble_names(style$name, times_gas)
      
      # names and style order for p-value figs
      myNames <- names(pval_tab)
      order <- match(myNames, style$name)
      
      for (myName in myNames) {
        # sort p-values, compare to uniform quantiles, calculate RMSD
        obj <- pvals_to_null_rmsd(pval_tab[[myName]], causal_indexes)
        # store RMSD value
        pval_rmsd[[myName]] <- obj$rmsd
        # find row in style, to get col/lty from
        pr <- pvals_to_pr(pval_tab[[myName]], causal_indexes)
        # store AUCs
        pr_auc[[myName]] <- pr$auc.integral
        
        
      }
      return(c( pval_rmsd, pr_auc))
    }
    
    
    #store mean of p and auc values after removing the casual index
    p<-pvals_auc_collect(name_out, causal_indexes, pvals, times_kin, times_gas, style=style)
    p<-unlist(p)
    M_auc[j,pcs]<-p[2]
    M_rmsd[j,pcs]<-p[1]
    
  }
  
  
}
setwd("/dscrhome/yy222/real_data_set")

write.table(M_auc, file="auc_k_fixed_k_10_pcs_2_30_n_100_10_3_human.txt")
write.table(M_rmsd, file="rmsd_k_fixed_k_10_pcs_2_30_n_100_10_3_human.txt")


