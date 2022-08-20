# this scripts runs PCA on simulation, with a certain number of PCs
# this version uses plink for the actual association, but uses my eigenvectors precomputed in R

library(optparse)
library(readr)
library(genbin) # binary wrappers

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
rep <- opt$rep
n_pcs <- opt$n_pcs
threads <- opt$threads
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit
env1 <- opt$env1
env2 <- opt$env2

# do this consistency check early
if ( !is.na( env1 ) && is.na( env2 ) )
    stop( 'If --env1 is specified, must also specify --env2!' )

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
# plink_glm will handle this correctly
if ( n_pcs == 0 )
    file_covar <- NULL

# adjust paths if using fes model or non-default m_causal_fac or herit
dir_phen <- '' # current dir
# use subdir instead in this case
if ( fes )
    dir_phen <- 'fes/'
if ( m_causal_fac != 10 )
    dir_phen <- paste0( dir_phen, 'm_causal_fac-', m_causal_fac, '/' )
if ( herit != 0.8 )
    dir_phen <- paste0( dir_phen, 'h', herit, '/' )
if ( !is.na( env1 ) )
    dir_phen <- paste0( dir_phen, 'env', env1, '-', env2, '/' )

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
data <- plink_glm(
    name_in_lower,
    name_phen = name_phen,
    name_out = name_out, # write outputs into current level, add number of PCs
    file_covar = file_covar,
    threads = threads
)

# save p-values into a file, simple human-readable format
write_lines( data$p, file_out )

# clean up when we're done with plink
# deletes GLM table and log only
# files are in current level (matches `name_out`)
delete_files_plink_glm( name_out )
delete_files_log( name_out )
