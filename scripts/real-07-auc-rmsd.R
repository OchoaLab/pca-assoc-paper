# this script scans for existing p-value calculations, and summarizes them into AUC and RMSD

library(optparse)
library(readr)
library(tibble)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('pvals_to_null_rmsd.R')
source('pvals_to_pr_auc.R')
setwd( dir_orig ) # go back to where we were

# constants
methods <- c('pca-plink', 'gcta')

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
                help = "Max replicates", metavar = "int")
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

for ( rep in 1 : rep_max ) {
    # move higher to the "reps" location
    # this is so GCTA's temporary files don't overwrite files from other parallel runs
    dir_out <- paste0( 'rep-', rep )
    setwd( dir_out )
    
    message( 'rep: ', rep )

    # load trait info we need
    load( 'simtrait.RData' )
    # loads:
    ## trait,
    ## causal_indexes,
    ## causal_coeffs

    # start a big loop
    for ( method in methods ) {
        message( 'method: ', method )
        for ( n_pcs in 1 : n_pcs_max ) {
            message( 'pcs: ', n_pcs )
            
            # file to read
            file_pvals <- paste0( 'pvals_', method, '_', n_pcs, '.txt.gz' )
            file_sum <- paste0( 'sum_', method, '_', n_pcs, '.txt.gz' )

            # if output is already there, don't do anything (don't recalculate)
            if ( file.exists( file_sum ) ) {
                message( 'File already exists, skipping: ', file_sum )
            } else if ( !file.exists( file_pvals ) ) {
                message( 'No p-vals, skipping: ', file_pvals )
            } else {
                # read the file
                pvals <- as.numeric( read_lines( file_pvals ) )
                # calculate RMSD_p
                rmsdp <- pvals_to_null_rmsd(pvals, causal_indexes)$rmsd
                # calculate AUC_PR
                aucpr <- pvals_to_pr_auc(pvals, causal_indexes)
                # put everything into a tibble, with all the info we want conveniently in place
                tib <- tibble(
                    method = method,
                    pc = n_pcs,
                    rep = rep,
                    rmsd = rmsdp,
                    auc = aucpr
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
