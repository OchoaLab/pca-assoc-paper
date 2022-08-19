# this script scans for existing p-value calculations, and summarizes them into AUC and RMSD

library(optparse)
library(readr)
library(tibble)
library(genio)
library(parallel)
library(simtrait) # pval_srmsd, pval_aucpr, pval_infl

# constants
methods <- c('pca-plink-pure', 'gcta') # 'pca-plink', 
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
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double")
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
    
    # if output is already there, don't do anything (don't recalculate)
    if ( file.exists( file_sum ) )
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
    
    # calculate RMSD_p
    rmsdp <- pval_srmsd( pvals, causal_indexes )
    # calculate AUC_PR
    aucpr <- pval_aucpr( pvals, causal_indexes )
    # calculate inflation factor from p-values (map back to Chi-Sq)
    lambda <- pval_infl( pvals )
    # calculate traditional type I error and calibrated power stats
    type_1_err <- pval_type_1_err( pvals, causal_indexes )
    power_calib <- pval_power_calib( pvals, causal_indexes )

    # put everything into a tibble, with all the info we want conveniently in place
    tib <- tibble(
        method = method,
        pc = n_pcs,
        rep = rep,
        rmsd = rmsdp,
        auc = aucpr,
        lambda = lambda,
        type_1_err = type_1_err,
        power_calib = power_calib
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
