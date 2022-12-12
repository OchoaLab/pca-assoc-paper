# this script scans for existing p-value calculations, and summarizes them into AUC and RMSD

library(optparse)
library(readr)
library(tibble)
library(genio)
library(parallel)
library(simtrait) # pval_srmsd, pval_aucpr, pval_infl

# constants
methods <- c('pca-plink-pure', 'gcta', 'gcta-labs') # 'pca-plink', 
# the name is for dir only, actual file is just "data"
name_in <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--n_pcs", type = "integer", default = 90,
                help = "Max number of PCs", metavar = "int"),
    make_option(c("-r", "--rep"), type = "integer", default = 50,
                help = "Max replicates", metavar = "int"),
    make_option("--sim", action = "store_true", default = FALSE, 
                help = "Genotypes are simulated (rather than real; alters location only)"),
    make_option(c("-t", "--threads"), type = "integer", default = detectCores(), 
                help = "number of threads (default use all cores, but may consume excessive memory)", metavar = "int"),
    make_option("--force", action = "store_true", default = FALSE, 
                help = "Overwrite all outputs (default is to keep existing outputs)"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option("--env1", type = "double", default = NA,
                help = "Variance of 1st (coarsest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option("--env2", type = "double", default = NA,
                help = "Variance of 2nd (finest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep_max <- opt$rep
n_pcs_max <- opt$n_pcs
threads <- opt$threads
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit
env1 <- opt$env1
env2 <- opt$env2

# do this consistency check early
if ( !is.na( env1 ) && is.na( env2 ) )
    stop( 'If --env1 is specified, must also specify --env2!' )

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

# in simulations, each rep has its own bim file, but they all have the same number of loci, so let's just read it from rep-1
if (opt$sim)
    name_in <- paste0( 'rep-1/', name_in )
# get number of loci
m_loci <- count_lines( name_in, 'bim' )

# function to process one n_pcs
# note every other parameter is a global variable
auc_rmsd_one_pcs <- function(n_pcs) {
    
    # file to read
    file_pvals <- paste0( 'pvals_', method, '_', n_pcs, '.txt.gz' )
    file_sum <- paste0( 'sum_', method, '_', n_pcs, '.txt.gz' )
    
    # if output is already there and not forcing overwrite, don't do anything (don't recalculate)
    if ( !opt$force && file.exists( file_sum ) )
        return()
    # if there's no input file, also skip silently
    if ( !file.exists( file_pvals ) )
        return()
    
    # only report files processed in this run
    message(
        'rep: ', rep,
        ', method: ', method,
        ', pcs: ', n_pcs
    )
    
    # read the file
    pvals <- read_lines(
        file_pvals,
        na = 'NA' # use this to prevent warnings
    )

    # special case is a single line with the word NULL, which occurs particularly with GCTA outputs when model does not converge/is underdetermined
    if ( length(pvals) == 1 && pvals == 'NULL' ) {
        pvals <- NULL
    } else {
        # make sure length is correct!
        # since we're in parallelized code, skip instead of dying (return early without generating output)
        if ( length(pvals) != m_loci ) {
            message( 'File has ', length(pvals), ' p-values, expected ', m_loci, '.  SKIPPING!' )
            return()
        }
        # convert strings to numeric
        pvals <- as.numeric( pvals )
    }

    # calculate type I error and calibrated power at a range of significance values
    alpha <- c(1e-2, 1e-4, 1e-6, 1e-8)
    type_1_err  <- pval_type_1_err(  pvals, causal_indexes, alpha = alpha )
    power_calib <- pval_power_calib( pvals, causal_indexes, alpha = alpha )
    
    # calculate all other desired statistics, put everything into a tibble
    tib <- tibble(
        method       = method,
        pc           = n_pcs,
        rep          = rep,
        rmsd         = pval_srmsd( pvals, causal_indexes ),
        auc          = pval_aucpr( pvals, causal_indexes ),
        lambda       = pval_infl( pvals ),
        type_1_err2  = type_1_err[1],
        type_1_err4  = type_1_err[2],
        type_1_err6  = type_1_err[3],
        type_1_err8  = type_1_err[4],
        power_calib2 = power_calib[1],
        power_calib4 = power_calib[2],
        power_calib6 = power_calib[3],
        power_calib8 = power_calib[4]
    )
    # save this one row to output file
    write_tsv(
        tib,
        file_sum
    )
}

# remember this location, using global path, for easily returning across multiple levels when done
dir_base <- getwd()

for ( rep in 1 : rep_max ) {
    # specify location of files to process, as many levels as needed
    dir_out <- paste0( 'rep-', rep, '/' )
    if ( fes )
        dir_out <- paste0( dir_out, 'fes/' )
    if ( m_causal_fac != 10 )
        dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
    if ( herit != 0.8 )
        dir_out <- paste0( dir_out, 'h', herit, '/' )
    if ( !is.na( env1 ) )
        dir_out <- paste0( dir_out, 'env', env1, '-', env2, '/' )
    
    # skip reps that we haven't calculated at all
    if ( !dir.exists( dir_out ) )
        next
    # else move to that destination
    setwd( dir_out )
    
    # load trait info we need
    load( 'simtrait.RData' )
    # loads:
    ## trait,          # not used
    ## causal_indexes, # only one used
    ## causal_coeffs   # not used

    # start a big loop
    for ( method in methods ) {
        # this is where parallelization happens!
        mclapply(
            0 : n_pcs_max,
            auc_rmsd_one_pcs,
            mc.cores = threads
        )
    }
    
    # return to base when done with this rep
    setwd( dir_base )
}
