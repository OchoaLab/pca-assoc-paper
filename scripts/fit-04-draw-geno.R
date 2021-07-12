# a real-data/tree version of sim-01-draw-geno.R

# creates:
# - data.bed
# - data.bim
# - data.fam
# - p_anc.RData

library(optparse)
library(bnpsd)
library(genio)
library(readr)

# hardcoded params
verbose <- TRUE # to print messages
# the name is for dir only, actual file is just "data"
name_out <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 100000, 
                help = "number of loci", metavar = "int"),
    make_option(c("-r", "--rep"), type = "integer", default = 1,
                help = "Replicate number", metavar = "int"),
    make_option("--maf_real", action = "store_true", default = FALSE, 
                help = "Draw ancestral allele frequencies from real distribution"),
    make_option("--maf_min", type = "double", default = 0.01, 
                help = "Minimum MAF, to simulate MAF-based ascertainment bias", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
m_loci <- opt$m_loci
rep <- opt$rep
maf_real <- opt$maf_real
maf_min <- opt$maf_min

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

##############################
### LOAD PRECALCUATED DATA ###
##############################

# directory above must already exist
setwd( '../data/' )

p_anc_distr <- NULL
if ( maf_real ) {
    # deduce original file path
    name_real <- sub( '_sim$', '', name )
    setwd( name_real )
    
    # load MAF
    load( 'maf.RData' )
    # loads: maf

    # go back down so we can find the simulated path
    setwd( '..' )
}

setwd( name )

# load precalculated bnpsd data from RData file
load( 'bnpsd.RData' )
# loads: admix_proportions, tree_subpops, fam

if (maf_real) {
    # need mean kinship of this simulation, which has been calculated already
    if ( name == 'HoPacAll_ld_prune_1000kb_0.3_maf-0.01_sim' ) {
        kinship_mean <- as.numeric( read_lines( 'kinship_mean.txt' ) )
    } else {
        # override with this hardcoded value that performed best in simulations
        kinship_mean <- 0.4
    }
    message( 'kinship_mean used for undiff_af: ', kinship_mean )
    
    # undifferentiate whole original distribution
    p_anc_distr <- undiff_af( maf, kinship_mean )$p
}

#####################
### SIM GENOTYPES ###
#####################

# create a new directory to specify replicate
dir_out <- paste0( 'rep-', rep )
# let's not overwrite things, under the assumption that the simulations are very expensive to redo
# if the output directory exists, assume all the files we want are there too.  Only a non-existent output directory will work
if ( dir.exists( dir_out ) ) {
    stop('Output exists, will not overwrite: ', dir_out)
} else {
    # create directory and move into there
    dir.create( dir_out )
    setwd( dir_out )
}

# draw allele freqs and genotypes
if (verbose)
    message('draw_all_admix')
# simulate data from this
# only difference is here number of loci will be smaller than real data, on purpose (and their MAF distribution will be practically uniform too)
out <- draw_all_admix(
    admix_proportions = admix_proportions,
    tree_subpops = tree_subpops,
    m_loci = m_loci,
    p_anc_distr = p_anc_distr,
    maf_min = maf_min,
    verbose = verbose
)
X <- out$X
p_anc <- out$p_anc

# write to plink BED/BIM/FAM
# use real FAM here, also ensures that X and fam agree in individual order
write_plink(
    name_out,
    X,
    fam = fam,
    verbose = verbose
)

# write final p_anc
save( p_anc, file = 'p_anc.RData' )
