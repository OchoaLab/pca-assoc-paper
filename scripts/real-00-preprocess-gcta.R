library(optparse)
library(genbin) # binary wrappers

# number of PCs to explore
n_pcs_max <- 90
# the name is for dir only, actual file is just "data"
name_in <- 'data'
# output for PCA holds number of PCs
name_out <- paste0( name_in, '-n_pcs_', n_pcs_max )

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character")
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

# main steps!
# GRM has same name as plink files
gcta_grm( name_in )
# eigenvec files have number of PCs suffixed to name
gcta_pca( name_in, name_out, n_pcs = n_pcs_max )

# cleanup
delete_files_log( name_in )
delete_files_log( name_out )
