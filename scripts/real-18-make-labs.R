# This script makes a categorical covariates file to be used for a narrow test, for GCTA.
# To test if direct use of true labels of env effects improves over PCs a lot or not.
# File is identical across replicates for each dataset (sim and real), so makes sense to just make once.
# Regression uses finest labels only (second level of simtrait versions), since coarse labels are special case of fine groupings in all cases.

library(optparse)
library(genio)
library(readr)

# the name is for dir only, actual file is just "data"
name_in <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--sim", action = "store_true", default = FALSE, 
                help = "Genotypes are simulated (rather than real; alters paths only)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

# simulations have fam in reps, just take first one (only thing that matters is number of individuals in this case, which is the same)
name_fam <- if ( opt$sim ) paste0( 'rep-1/', name_in ) else name_in

# load FAM table, for phen output
fam <- read_fam( name_fam )

# define `labs`!
if ( opt$sim ) {
    # there's no labels here, but use geographic order for grouping individuals
    k_subpops <- 25 # hardcoded!  n/k =  40
    n_ind <- nrow( fam )
    labs <- ceiling( ( 1 : n_ind ) * k_subpops / n_ind )
} else {
    # these are already aligned in fam file!
    labs <- fam$fam
}

# start from fam to write output
# just add labels as a new column
fam$lab <- labs

# write this!
# uses a genio internal to have its nice features (subsets and reorders columns, performs checks, writes no header as required by GCTA, etc)
genio:::write_tab_generic(
    file = name_in,
    tib = fam,
    ext = 'covar',
    tib_names = c('fam', 'id', 'lab')
)
