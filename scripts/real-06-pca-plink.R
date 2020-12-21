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
    make_option("--plink", action = "store_true", default = FALSE, 
                help = "use PCs calculated with plink (default is to use PCs from `kinship_std`)"),
    make_option("--sim", action = "store_true", default = FALSE, 
                help = "Genotypes are simulated (rather than real; alters paths only)"),
    make_option("--dcc", action = "store_true", default = FALSE, 
                help = "Duke Compute Cluster runs (alters paths only)"),
    make_option(c("-t", "--threads"), type = "integer", default = 0, 
                help = "number of threads (default use all cores if not DCC, 1 core if DCC)", metavar = "int"),
    make_option("--const_herit_loci", action = "store_true", default = FALSE, 
                help = "Causal coefficients constructed to result in constant per-locus heritability (saved in diff path)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep <- opt$rep
n_pcs <- opt$n_pcs
threads <- opt$threads
const_herit_loci <- opt$const_herit_loci

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

# move higher to the "reps" location
# this is so GCTA's temporary files don't overwrite files from other parallel runs
dir_out <- paste0( 'rep-', rep )
setwd( dir_out )

###########
### PCA ###
###########

# identify the PCs used
# (std are from kinship_std, plink are for pure plink like others would)
# for input covariates file and for a temporary output path
name_pcs <- if ( opt$plink ) 'plink' else 'std'

# this is for final outputs
method <- if ( opt$plink ) 'pca-plink-pure' else 'pca-plink'

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

file_covar <- paste0( name_in_lower, '-', name_pcs, '-n_pcs_', n_pcs, '.eigenvec' )
# there's no covariates file to pass if we want zero PCs
# gas_plink will handle this correctly
if ( n_pcs == 0 )
    file_covar <- NULL

# adjust paths if using const_herit_loci model
dir_phen <- '' # current dir
# use subdir instead in this case
if ( const_herit_loci )
    dir_phen <- 'const_herit_loci/'

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
    threads = threads,
    ver_older = opt$dcc # prevents an error when there are no covariates
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
