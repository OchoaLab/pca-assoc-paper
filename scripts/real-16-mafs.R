# gathers MAF distributions, for a plot of interest

library(BEDMatrix)
library(simtrait) # for allele_freqs
library(optparse)

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
                help = "Genotypes are simulated (rather than real; alters paths only, reads and writes data in 'rep-1/' subdir)")
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

# if data is simulated, genotypes are not in base dir but in rep-*/ subdirs (one for every rep)
# all reps ought to be similar, just use rep-1 (no need to ever have to repeat this)
if ( opt$sim )
    setwd( 'rep-1' )

# load genotypes BEDMatrix object
X <- BEDMatrix( name_in )

# and get back down when done
if ( opt$sim )
    setwd( '..' )

# calculate desired allele frequencies
maf <- allele_freqs( X, fold = TRUE )

# save all data on `data` base directory
save( maf, file = 'maf.RData' )
