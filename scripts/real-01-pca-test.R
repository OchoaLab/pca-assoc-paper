library(optparse)
library(popkinsuppl) # for PCA's kinship estimator
library(BEDMatrix)
library(readr)
library(dplyr)
library(popkin)
library(ochoalabtools)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('kinship_to_evd.R')
setwd( dir_orig ) # go back to where we were

# number of PCs to explore
n_pcs_max <- 90
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

# move to where the data is
setwd( '../data/' )
setwd( name )

#################
### LOAD GCTA ###
#################

# first load GCTA's precomputed PCs, for comparison
pcs_gcta <- read_table2(
    paste0( name_in, '-n_pcs_', n_pcs_max, '.eigenvec' ),
    col_names = FALSE
)
# first two columns (FAM/ID) are not needed
# but copy them so we can write table in same format
fam <- pcs_gcta[ , 1:2 ]
pcs_gcta <- pcs_gcta[ , -(1:2) ]
# turn to numeric matrix
pcs_gcta <- as.matrix( pcs_gcta )

##############
### MY PCA ###
##############

# load with BEDMatrix
X <- BEDMatrix( name_in )

# estimate kinship (std) and eigenvectors for PCA
# estimate kinship the old way
kinship_estimate_old <- kinship_std(X)
# get all eigenvalues
eigenvectors_estimate_old <- kinship_to_evd( kinship_estimate_old )
# subset to match GCTA's (and all we need for our tests, in any case)
eigenvectors_estimate_old <- eigenvectors_estimate_old[ , 1 : n_pcs_max ]

# add columns that match GCTA's, to try things both ways
tib <- bind_cols( fam, as_tibble( eigenvectors_estimate_old ) )
# save in a file like GCTA's
write_tsv(
    tib,
    paste0( name_in, '-std-n_pcs_', n_pcs_max, '.eigenvec' ),
    col_names = FALSE
)

##################
### COMPARISON ###
##################

# the key comparison is this inner product of vectors, which in a perfect agreement leads to the identity matrix (as these are orthonormal)
# results in a dim-`n_pcs_max` square matrix (i.e. 90 x 90)
comparison <- crossprod( pcs_gcta, eigenvectors_estimate_old )

# a crude comparison shows that only the top eigenvectors sort of agree
# this is ok as r becomes very large, they don't fit well regardless, but some of the reranking is surprising
fig_start(
    'pca-comparison',
    width = 4,
    mar_b = 0,
    mar_l = 0
)
plot_popkin(
    abs( comparison ),
    ylab = 'Eigenvector',
    leg_title = '| Inner Product |'
)
fig_end()
