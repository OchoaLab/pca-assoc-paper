# simulates and writes genotypes to plink BED/BIM/FAM files
# based closely on ~/scripts/sim_reps_01_sim_geno.R on main scripts dir of this repo

# creates:
# - data.bed
# - data.bim
# - data.fam
# - p_anc.RData

library(optparse) # for terminal options
library(bnpsd)    # simulate admixed population (structured population)
library(genio)    # to write BED files for external software
library(tibble)   # for initializing fake bim

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
# sim_children_generations_* code
source('draw_geno_child.R')
setwd( dir_orig ) # go back to where we were

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
generations <- opt$generations
rep <- opt$rep

# name of directory containing data
name <- paste0(
    'sim',
    '-n', n_ind,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-g', generations
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

# load parents (pedigree) info if needed
if (generations > 1) {
    load( 'parents.RData' )
    # loads: parents
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
out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci)
X <- out$X # genotypes
p_anc <- out$p_anc # ancestral AFs

if (generations > 1) {
    # simulate realistic generations of families

    if (verbose)
        message('sim_children_generations_genotypes')
    # final genotypes
    X <- sim_children_generations_genotypes(X, parents, verbose = verbose)

    # NOTE: p_anc doesn't change, this is accurate
    
    # handle fixed loci (a big pain!)
    fixed_loci_indexes <- fixed_loci(X)
    m_loci_fixed <- sum( fixed_loci_indexes )
    while (m_loci_fixed > 0) { # draw things anew, over and over until nothing was fixed
        # draw allele freqs and genotypes
        out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci_fixed)
        # overwrite fixed loci with redrawn polymorphic data
        #            X[fixed_loci_indexes, ] <- out$X # genotypes
        p_anc[fixed_loci_indexes] <- out$p_anc # ancestral AFs

        if (verbose)
            message('sim_children_generations_genotypes (redrawn)')
        # repeat children draws through generations, then
        # overwrite fixed loci with redrawn (hopefully) polymorphic data
        X[fixed_loci_indexes, ] <- sim_children_generations_genotypes(out$X, parents, verbose = verbose)
        
        # look for remaining fixed loci (to continue or stop loop)
        fixed_loci_indexes <- fixed_loci(X)
        m_loci_fixed <- sum( fixed_loci_indexes )
    }
}

# write to plink BED/BIM/FAM 
write_plink(
    name_out,
    X
)

# write final p_anc
save( p_anc, file = 'p_anc.RData' )
