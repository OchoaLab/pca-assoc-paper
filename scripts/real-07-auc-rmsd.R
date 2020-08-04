# this script scans for existing p-value calculations, and summarizes them into AUC and RMSD

library(optparse)
library(readr)
library(tibble)
library(genio)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('pvals_to_null_rmsd.R')
source('pvals_to_pr_auc.R')
setwd( dir_orig ) # go back to where we were

# constants
methods <- c('pca-plink', 'gcta')
# the name is for dir only, actual file is just "data"
name_in <- 'data'
# constant factor needed to transform median p-values into inflation factors lambda from chi-square
df <- 1
x_m <- qchisq( 0.5, df = df )

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
                help = "Genotypes are simulated (rather than real; alters location only)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep_max <- opt$rep
n_pcs_max <- opt$n_pcs

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

for ( rep in 1 : rep_max ) {
    # move higher to the "reps" location
    # this is so GCTA's temporary files don't overwrite files from other parallel runs
    dir_out <- paste0( 'rep-', rep )
    # skip reps that we haven't calculated at all
    if ( !dir.exists( dir_out ) )
        next
    setwd( dir_out )
    
    # load trait info we need
    load( 'simtrait.RData' )
    # loads:
    ## trait,
    ## causal_indexes,
    ## causal_coeffs

    # start a big loop
    for ( method in methods ) {
        for ( n_pcs in 0 : n_pcs_max ) {
            # file to read
            file_pvals <- paste0( 'pvals_', method, '_', n_pcs, '.txt.gz' )
            file_sum <- paste0( 'sum_', method, '_', n_pcs, '.txt.gz' )

            # if output is already there, don't do anything (don't recalculate)
            if ( file.exists( file_sum ) ) {
                ## message( 'File already exists, skipping: ', file_sum )
            } else if ( !file.exists( file_pvals ) ) {
                ## message( 'No p-vals, skipping: ', file_pvals )
            } else {
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
                    # when pvals are NULL, they get printed as an empty file
                    # this is what we get from pvals_to_null_rmsd, pvals_to_pr_auc if we had set pvals==NULL
                    rmsdp <- NA
                    aucpr <- NA
                    # this special case wasn't handled below, I think
                    lambda <- NA
                } else {
                    # make sure length is correct!
                    if ( length(pvals) != m_loci )
                        stop( 'File has ', length(pvals), ' p-values, expected ', m_loci )

                    # convert strings to numeric
                    pvals <- as.numeric( pvals )
                    # calculate RMSD_p
                    rmsdp <- pvals_to_null_rmsd(pvals, causal_indexes)$rmsd
                    # calculate AUC_PR
                    aucpr <- pvals_to_pr_auc(pvals, causal_indexes)
                    # calculate inflation (but on p-values instead of Chi-Sq; this makes sense as not all stats are Chi-Sq anyway)
                    # NOTES:
                    # - not subsetting for true nulls (as done in practice)
                    # - exact inversion to chi-sq stat
                    lambda <- qchisq( 1 - median(pvals, na.rm = TRUE), df = df ) / x_m
                }
                # put everything into a tibble, with all the info we want conveniently in place
                tib <- tibble(
                    method = method,
                    pc = n_pcs,
                    rep = rep,
                    rmsd = rmsdp,
                    auc = aucpr,
                    lambda = lambda
                )
                # save this one row to output file
                write_tsv(
                    tib,
                    file_sum
                )
            }
        }
    }
    
    # move back down when done with this rep
    setwd( '..' )
}
