# compute PCs using plink, figure out how they relate to PCs from my method (kinship_std) and GCTA's

library(optparse)
library(popkin) # for comparison plot only
library(genio) # for read_eigenvec
library(ochoalabtools)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('gas_plink.R')
source('paths.R')
setwd( dir_orig ) # go back to where we were

# number of PCs to explore
n_pcs_max <- 90
# the name is for dir only, actual file is just "data"
name_in <- 'data'
# plink: here there's no timing, so let's be as fast as possible (use all threads)
threads <- 0
# plink need a MAF threshold, otherwise their standard kinship (MOR version) behaves poorly
maf <- 0.1

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

#############
### PLINK ###
#############

# this is presumably the approach most people use

# unique name for this run
name_out <- paste0( name_in, '-plink-n_pcs_', n_pcs_max )

# this creates the PCA file
gas_plink_pca(
    plink2_bin,
    name = name_in,
    name_out = name_out,
    n_pcs = n_pcs_max,
    maf = maf,
    threads = threads
)

# cleanup
invisible( file.remove( paste0( name_out, '.log' ) ) )

# load the data for comparisons
out <- read_eigenvec(
    file = name_out
)
eigenvec_plink <- out$eigenvec

#################
### LOAD GCTA ###
#################

# first load GCTA's precomputed PCs, for comparison
out <- read_eigenvec(
    file = paste0( name_in, '-n_pcs_', n_pcs_max )
)
eigenvec_gcta <- out$eigenvec

## ###################
## ### LOAD MY PCA ###
## ###################

## # also load my custom precomputed PCs (from kinship_std), for comparison
## out <- read_eigenvec(
##     file = paste0( name_in, '-std-n_pcs_', n_pcs_max )
## )
## eigenvec_std <- out$eigenvec

###################
### COMPARISONS ###
###################

# the key comparison is this inner product of vectors, which in a perfect agreement leads to the identity matrix (as these are orthonormal)
# results in a dim-`n_pcs_max` square matrix (i.e. 90 x 90)
# in this version we do 3 times!
#comparison_gcta_std <- crossprod( eigenvec_gcta, eigenvec_std )
comparison_gcta_plink <- crossprod( eigenvec_gcta, eigenvec_plink )
#comparison_std_plink <- crossprod( eigenvec_std, eigenvec_plink )

# a crude comparison shows that only the top eigenvectors sort of agree
# this is ok as r becomes very large, they don't fit well regardless, but some of the reranking is surprising
fig_start(
    'pca-comparison2',
    width = 4, # 10,
    mar_t = 2,
    mar_b = 0,
    mar_l = 0
)
plot_popkin(
    list(
#        abs( comparison_gcta_std ),
        abs( comparison_gcta_plink )
#        abs( comparison_std_plink )
    ),
    titles = c(
#        'GCTA-STD',
        'GCTA-plink'
#        'STD-plink'
    ),
    ylab = 'Eigenvector',
    leg_title = 'Abs. Correlation'
)
fig_end()
