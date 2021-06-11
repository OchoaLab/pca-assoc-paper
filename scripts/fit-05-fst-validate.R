library(optparse)
library(popkin)
library(genio)
library(readr)
library(ochoalabtools)
library(tibble)

# estimates FST from actual genotype simulated data, for validation

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

# need original annotations, get from source real data file
# assuming a pattern here:
name_real <- sub( '_sim$', '', name )

# move to where the data is
setwd( '../data/' )

# first load annotations
setwd( name_real )
subpop_info <- read_tsv('pops-annot.txt', comment = '#')

# now to where the simulated data is
setwd( '..' )
setwd( name )
# for this data, fam and kinship are in rep-1 only
# will save FST here too only
setwd( 'rep-1' )

# read full fam data
fam <- read_fam( name_in )

# read existing popkin data
data <- read_grm( 'popkin' )
kinship <- data$kinship

# map subpopulations using sub-subpopulations
fam$superpop <- subpop_info$superpop[ match( fam$fam, subpop_info$pop ) ]

# confirm matching orders here too
stopifnot( fam$id == colnames( kinship ) )
## # reorder individuals so subpopulations come out in desired order:
## indexes <- order( match( fam$fam, subpop_info$pop ) )
## fam <- fam[ indexes, ]
## kinship <- kinship[ indexes, indexes ]

# calculate weights that balance everything!
weights <- weights_subpops( fam$superpop, fam$fam )
fst_est <- fst( kinship, weights )
# save FST for table
write_lines( fst_est, 'fst_est.txt' )
