# simulates and writes genotypes to plink BED/BIM/FAM files
# based closely on ~/scripts/sim_reps_01_sim_geno.R on main scripts dir of this repo

# creates:
# - data.bed
# - data.bim
# - data.fam
# - p_anc.RData
# - all.fam (if generations > 1)

library(optparse) # for terminal options
library(bnpsd)    # simulate admixed population (structured population)
library(genio)    # to write BED files for external software
library(tibble)   # for initializing fake bim
library(simfam)   # to simulate family structure

# hardcoded params
verbose <- TRUE # to print messages
# the name is for dir only, actual file is just "data"
name_out <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-n", "--n_ind"), type = "integer", default = 1000,
                help = "number of individuals", metavar = "int"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 100000, 
                help = "number of loci", metavar = "int"),
    make_option(c("-k", "--k_subpops"), type = "integer", default = 10, 
                help = "admixture intermediate subpopulations", metavar = "int"),
    make_option(c("-f", "--fst"), type = "double", default = 0.1, 
                help = "FST (fixation index)", metavar = "double"),
    make_option(c("--bias_coeff"), type = "double", default = 0.5, 
                help = "admixture bias coeff", metavar = "double"),
    make_option(c("-g", "--generations"), type = "integer", default = 1, 
                help = "number of generations, for realistic local kinship", metavar = "int"),
    make_option(c("-r", "--rep"), type = "integer", default = 1,
                help = "Replicate number", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_ind <- opt$n_ind
m_loci <- opt$m_loci
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
G <- opt$generations
rep <- opt$rep

# name of directory containing data
name <- paste0(
    'sim',
    '-n', n_ind,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-g', G
)

# directory above must already exist
setwd( '../data/' )
setwd( name )

##############################
### LOAD PRECALCUATED DATA ###
##############################

# load precalculated bnpsd data from RData file
load( 'bnpsd.RData' )
# loads: admix_proportions, inbr_subpops

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
out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci)
X <- out$X # genotypes
p_anc <- out$p_anc # ancestral AFs

# leave undefined if there's no family simulation
fam <- NULL

if (G > 1) {
    # simulate realistic generations of families

    # first simulate pedigree (varies per replicate)
    if (verbose)
        message('sim_pedigree')
    data_simfam <- sim_pedigree( n_ind, G )
    fam <- data_simfam$fam
    ids <- data_simfam$ids
    # save as ordinary FAM file!
    # NOTE: this contains individuals without descendants and also all previous generations,
    # so it's not redundant with the final output `data.fam` (only has last generation)
    write_fam( 'all.fam', fam, verbose = verbose )

    # prune fam table now, to not simulate unnecessary individuals without descendants
    fam <- prune_fam( fam, ids[[ G ]] )

    # `simfam` requires names for `X`
    colnames( X ) <- ids[[ 1 ]]
    
    if (verbose)
        message('geno_last_gen')
    # replace genotypes of founders with genotypes for last generation
    X <- geno_last_gen( X, fam, ids )
    
    # NOTE: p_anc doesn't change, this is accurate
    
    # handle fixed loci (a big pain!)
    fixed_loci_indexes <- fixed_loci(X)
    m_loci_fixed <- sum( fixed_loci_indexes )
    while (m_loci_fixed > 0) { # draw things anew, over and over until nothing was fixed
        # draw allele freqs and genotypes
        out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci_fixed)
        # overwrite fixed loci with redrawn polymorphic data
        p_anc[fixed_loci_indexes] <- out$p_anc # ancestral AFs
        X_redrawn <- out$X # renamed for clarity
        # `simfam` requires names for `X_redrawn`
        colnames( X_redrawn ) <- ids[[ 1 ]]

        if (verbose)
            message('geno_last_gen (redrawn)')
        # repeat children draws through generations, then
        # overwrite fixed loci with redrawn (hopefully) polymorphic data
        X[fixed_loci_indexes, ] <- geno_last_gen( X_redrawn, fam, ids )
        
        # look for remaining fixed loci (to continue or stop loop)
        fixed_loci_indexes <- fixed_loci(X)
        m_loci_fixed <- sum( fixed_loci_indexes )
    }

    # replace `fam` (all generations) with last generation data only
    fam <- fam[ fam$id %in% ids[[G]], ]
}

# write to plink BED/BIM/FAM 
write_plink(
    name_out,
    X,
    fam = fam
)

# write final p_anc
save( p_anc, file = 'p_anc.RData' )
