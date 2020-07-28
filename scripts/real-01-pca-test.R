library(optparse)
library(popkinsuppl) # for PCA's kinship estimator
library(BEDMatrix)
library(popkin)
library(genio)
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

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

#################
### LOAD GCTA ###
#################

# first load GCTA's precomputed PCs, for comparison
out <- read_eigenvec(
    file = paste0( name_in, '-n_pcs_', n_pcs_max )
)
fam <- out$fam
eigenvec_gcta <- out$eigenvec

##############
### MY PCA ###
##############

# load with BEDMatrix
X <- BEDMatrix( name_in )

# estimate kinship (std) and eigenvectors for PCA
# estimate kinship the old way
kinship_estimate_old <- kinship_std(X)
# get all eigenvalues
eigenvec_std <- kinship_to_evd( kinship_estimate_old )
# subset to match GCTA's (and all we need for our tests, in any case)
eigenvec_std <- eigenvec_std[ , 1 : n_pcs_max ]

# save in a file like GCTA's
write_eigenvec(
    file = paste0( name_in, '-std-n_pcs_', n_pcs_max ),
    eigenvec = eigenvec_std,
    fam = fam
)

##################
### COMPARISON ###
##################

# the key comparison is this inner product of vectors, which in a perfect agreement leads to the identity matrix (as these are orthonormal)
# results in a dim-`n_pcs_max` square matrix (i.e. 90 x 90)
comparison <- crossprod( eigenvec_gcta, eigenvec_std )

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
