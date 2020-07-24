library(optparse)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('gas_lmm_gcta.R')
source('paths.R')
setwd( dir_orig ) # go back to where we were

# number of PCs to explore
n_pcs_max <- 90
# not timing this, so use all threads
threads <- 0
# the name is for dir only, actual file is just "data"
name_in <- 'data'

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

message("GCTA (GRM)")
gas_lmm_gcta_kin(gcta_bin, name_in, threads = threads)

message("GCTA (PCA)")
gas_lmm_gcta_pca(gcta_bin, name_in, threads = threads, n_pcs = n_pcs_max)
