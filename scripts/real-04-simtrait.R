# this scripts simulates a random trait for the given real dataset, storing key values

library(optparse)
library(simtrait)
library(genio)
library(BEDMatrix)

# the name is for dir only, actual file is just "data"
name_in <- 'data'
# minimum MAF for causal loci
min_maf_causal <- 0.01

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option(c("-r", "--rep"), type = "integer", default = 1,
                help = "Replicate number", metavar = "int"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--sim", action = "store_true", default = FALSE, 
                help = "Genotypes are simulated (rather than real; alters paths only)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
herit <- opt$herit
m_causal_fac <- opt$m_causal_fac
rep <- opt$rep
fes <- opt$fes

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )
dir_out <- paste0( 'rep-', rep )
if ( opt$sim )
    setwd( dir_out )

################
### simtrait ###
################

# set p_anc and kinship, as needed for sim and real cases
if ( opt$sim ) {
    # load true ancestral allele frequencies of simulation
    load( 'p_anc.RData' )
    # loads: p_anc
    kinship <- NULL
} else {
    # load precalculated popkin kinship matrix
    kinship <- read_grm( 'popkin' )$kinship
    p_anc <- NULL
}

# load FAM table, for phen output
fam <- read_fam( name_in )

# now that we have number of individuals, use a tenth of that for m_causal, rounding
m_causal <- round( nrow( fam ) / m_causal_fac )
message( 'm_causal: ', m_causal )

# load genotype data with BEDMatrix
X <- BEDMatrix( name_in, simple_names = TRUE )

# simulate trait
# sim genotypes use p_anc version, real data use kinship version
obj_trait <- sim_trait(
    X = X,
    m_causal = m_causal,
    herit = herit,
    p_anc = p_anc,
    kinship = kinship,
    fes = fes,
    maf_cut = min_maf_causal
)
# extract data of interest
# trait vector
trait = obj_trait$trait
# randomly-picked causal locus index
causal_indexes = obj_trait$causal_indexes
# locus coefficient vector (for causal loci only)
causal_coeffs = obj_trait$causal_coeffs

#############
### WRITE ###
#############

# in simulation, rep-*/ already exist and we're already in it
# in real data, rep-*/ doesn't exist, so create then move in there
if ( !opt$sim ) {
    if( !dir.exists( dir_out ) )
        dir.create( dir_out )
    setwd( dir_out )
}

# specify location of outputs, as many levels as needed
dir_out <- ''
if ( fes )
    dir_out <- paste0( dir_out, 'fes/' )
if ( m_causal_fac != 10 )
    dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
if ( herit != 0.8 )
    dir_out <- paste0( dir_out, 'h', herit, '/' )
# create if needed
if ( !dir.exists( dir_out ) )
    dir.create( dir_out, recursive = TRUE )
setwd( dir_out )

# check that outputs don't exist already
if ( file.exists( 'simtrait.RData' ) )
    stop( 'Output exists, will not overwrite: simtrait.RData' )
if ( file.exists( 'data.fam' ) )
    stop( 'Output exists, will not overwrite: data.fam' )

# now save, as R data
save(
    trait,
    causal_indexes,
    causal_coeffs,
    file = 'simtrait.RData'
)

# save trait as phen file too, for GCTA and plink
fam$pheno <- trait
write_phen( name_in, fam )
