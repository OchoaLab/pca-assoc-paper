# this scripts runs PCA on simulation, with a certain number of PCs
# this version uses plink for the actual association, but uses my eigenvectors precomputed in R

library(optparse)
library(readr)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('gas_plink.R')
source('paths.R')
setwd( dir_orig ) # go back to where we were

# the name is for dir only, actual file is just "data"
name_in <- 'data'
# identify the PCs used, for input covariates file and for a temporary output path
name_pcs <- 'plink'
# for final outputs
method <- 'pca-plink-pure'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--n_pcs", type = "integer", default = 0,
                help = "Number of PCs to use", metavar = "int"),
    make_option(c("-r", "--rep"), type = "integer", default = 1,
                help = "Replicate number", metavar = "int"),
    make_option("--sim", action = "store_true", default = FALSE, 
                help = "Genotypes are simulated (rather than real; alters paths only)"),
    make_option("--dcc", action = "store_true", default = FALSE, 
                help = "Duke Compute Cluster runs (alters paths only)"),
    make_option(c("-t", "--threads"), type = "integer", default = 0, 
                help = "number of threads (default use all cores if not DCC, 1 core if DCC)", metavar = "int"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep <- opt$rep
n_pcs <- opt$n_pcs
threads <- opt$threads
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac

# figure out what threads should be
# above default should be modified if --dcc and if the threads aren't zero (default, means use all, which we should never do on a cluster!)
if ( opt$dcc && threads == 0 )
    threads <- 1

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
if ( opt$dcc ) {
    # on DCC we go here, no name for simplicity (data is only temporarily there, so I won't have more than one dataset there at the time)
    setwd( '/work/ao128/' )
} else {
    setwd( '../data/' )
}
setwd( name )

# new level to this hierarchy
if ( m_causal_fac != 10 ) {
    dir_out <- paste0( 'm_causal_fac-', m_causal_fac )
    # now move in there
    setwd( dir_out )
}

# move higher to the "reps" location
# this is so GCTA's temporary files don't overwrite files from other parallel runs
dir_out <- paste0( 'rep-', rep )
setwd( dir_out )

###########
### PCA ###
###########

# message so we know where we're at
message(
    'rep: ', rep,
    ', method: ', method,
    ', pcs: ', n_pcs
)

# genotypes, PCs:
# - in real data, are all in lower level (shared across reps)
# - in simulated data, are all in current level (not shared across reps)
name_in_lower <- if ( opt$sim ) name_in else paste0( '../', name_in )
# data is even lower still in this mode
if ( !opt$sim && m_causal_fac != 10 )
    name_in_lower <- paste0( '../', name_in_lower )

file_covar <- paste0( name_in_lower, '-', name_pcs, '-n_pcs_', n_pcs, '.eigenvec' )
# there's no covariates file to pass if we want zero PCs
# gas_plink will handle this correctly
if ( n_pcs == 0 )
    file_covar <- NULL

# adjust paths if using fes model
dir_phen <- '' # current dir
# use subdir instead in this case
if ( fes )
    dir_phen <- 'fes/'

# only these are in dir_phen
name_phen <- paste0( dir_phen, name_in )
# output name for plink runs should have number of PCs, so concurrent runs don't overwrite each other
# not used in any final outputs (so only goal is to avoid concurrent run overlaps)
name_out <- paste0( dir_phen, 'plink_', name_pcs, '_', n_pcs )
# file to create
file_out <- paste0( dir_phen, 'pvals_', method, '_', n_pcs, '.txt.gz' )

# do not redo run if output was already present!
if ( file.exists( file_out ) )
    stop( 'Output already exists, skipping: ', file_out )

# actual run
obj <- gas_plink(
    plink2_bin,
    name = name_in_lower,
    name_phen = name_phen,
    name_out = name_out, # write outputs into current level, add number of PCs
    file_covar = file_covar,
    threads = threads
)

# all we care to preserve here are the p-values
pvals <- obj$pvals
# save into a file, simple human-readable format
write_lines(
    pvals,
    file_out
)

# clean up when we're done with plink
# deletes GAS table and log only
# files are in current level (matches earlier `name_out` option)
delete_files_plink_gas( name_out )
