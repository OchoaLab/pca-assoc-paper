# compute PCs using plink, figure out how they relate to PCs from GCTA

library(optparse)
library(genbin) # binary wrappers
#library(popkin) # for comparison plot only
#library(genio) # for read_eigenvec
#library(ochoalabtools)

# number of PCs to explore
n_pcs_max <- 90
# the name is for dir only, actual file is just "data"
name_in <- 'data'
# output for PCA holds number of PCs
name_out <- paste0( name_in, '-plink-n_pcs_', n_pcs_max )

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

# this creates the PCA file
plink_pca( name_in, name_out, n_pcs = n_pcs_max )

# cleanup
delete_files_log( name_out )

## # load the data for comparisons
## out <- read_eigenvec(
##     file = name_out
## )
## eigenvec_plink <- out$eigenvec

## #################
## ### LOAD GCTA ###
## #################

## # first load GCTA's precomputed PCs, for comparison
## out <- read_eigenvec(
##     file = paste0( name_in, '-n_pcs_', n_pcs_max )
## )
## eigenvec_gcta <- out$eigenvec

## ###################
## ### COMPARISONS ###
## ###################

## # the key comparison is this inner product of vectors, which in a perfect agreement leads to the identity matrix (as these are orthonormal)
## # results in a dim-`n_pcs_max` square matrix (i.e. 90 x 90)
## comparison_gcta_plink <- crossprod( eigenvec_gcta, eigenvec_plink )

## # a crude comparison shows that only the top eigenvectors sort of agree
## # this is ok as r becomes very large, they don't fit well regardless, but some of the reranking is surprising
## fig_start(
##     'pca-comparison2',
##     width = 4,
##     mar_t = 2,
##     mar_b = 0,
##     mar_l = 0
## )
## plot_popkin(
##     abs( comparison_gcta_plink ),
##     titles = 'GCTA-plink',
##     ylab = 'Eigenvector',
##     leg_title = 'Abs. Correlation'
## )
## fig_end()
