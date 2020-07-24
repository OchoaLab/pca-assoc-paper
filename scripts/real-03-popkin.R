# this scripts estimates kinship with popkin
# needed for mean kinship, to simulate trait heritability without bias

library(optparse)
library(popkin)
library(genio)
library(BEDMatrix)

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

##############
### popkin ###
##############

# load genotype data with BEDMatrix
X <- BEDMatrix( name_in, simple_names = TRUE )
# load FAM table, for subpopulation labels
fam <- read_fam( name_in )

# now estimate kinship
kinship <- popkin( X, fam$fam )

# save matrix for later use
write_grm( 'popkin', kinship, fam = fam )
