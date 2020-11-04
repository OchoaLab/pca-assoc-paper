# this scripts simulates a random trait for the given real dataset, storing key values
# NOTE: m_causal = n_ind / 10 is hardcoded

library(optparse)
library(simtrait)
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
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    ## make_option("--m_causal", type = "integer", default = 1000,
    ##             help = "num causal loci", metavar = "int"),
    make_option(c("-r", "--rep"), type = "integer", default = 1,
                help = "Replicate number", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
herit <- opt$herit
## m_causal <- opt$m_causal
rep <- opt$rep

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

################
### simtrait ###
################

# load precalculated popkin kinship matrix
kinship <- read_grm( 'popkin' )$kinship

# load FAM table, for phen output
fam <- read_fam( name_in )

# now that we have number of individuals, use a tenth of that for m_causal, rounding
m_causal <- round( nrow( fam ) / 10 )
message( 'm_causal: ', m_causal )

# load genotype data with BEDMatrix
X <- BEDMatrix( name_in, simple_names = TRUE )

# simulate trait
obj_trait <- sim_trait(
    X = X,
    m_causal = m_causal,
    herit = herit,
    kinship = kinship
)
# extract data of interest
# trait vector
trait = obj_trait$trait
# randomly-picked causal locus index
causal_indexes = obj_trait$causal_indexes
# locus effect size vector (for causal loci only)
causal_coeffs = obj_trait$causal_coeffs

#########################################
### write data in rep-specific subdir ###
#########################################

dir_out <- paste0( 'rep-', rep )
# create if it didn't already exist
if( !dir.exists( dir_out ) )
    dir.create( dir_out )
# now move in there
setwd( dir_out )

# now save, as R data
save(
    trait,
    causal_indexes,
    causal_coeffs,
    file = 'simtrait.RData'
)

# save trait as phen file too, for GCTA
fam$pheno <- trait
write_phen( name_in, fam )
