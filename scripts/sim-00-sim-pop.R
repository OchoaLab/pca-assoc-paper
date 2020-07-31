# simulates population structure (without genotypes)
# based closely on ~/scripts/sim_reps_00_sim_pop.R on main scripts dir of this repo

# things that can't be changed:
# - Fst of intermediate subpopulations is a ramp proportional to 1:k
# - Admixture proportions are from 1D geography model

# creates:
# - bnpsd.RData
# - parents.RData (if generations > 1)

library(optparse) # for terminal options
library(bnpsd)    # simulate admixed population (structured population)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
# sim_children_generations_* code
source('draw_geno_child.R')
setwd( dir_orig ) # go back to where we were

# hardcoded params
iterations <- 100 # GENERATIONS number of iterations until "pairing" code gives up and restarts pairing simulation at an earlier generation (in small populations there may be no solution as there are not enough sufficiently unrelated individuals)
verbose <- TRUE # to print messages
# the name is for dir only, actual file is just "data"
name_out <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-n", "--n_ind"), type = "integer", default = 1000, # CHANGE: or 1000
                help = "number of individuals", metavar = "int"),
    make_option(c("-k", "--k_subpops"), type = "integer", default = 10, 
                help = "admixture intermediate subpopulations", metavar = "int"),
    make_option(c("-f", "--fst"), type = "double", default = 0.1, 
                help = "FST (fixation index)", metavar = "double"),
    make_option(c("--bias_coeff"), type = "double", default = 0.5, 
                help = "admixture bias coeff", metavar = "double"),
    make_option(c("-g", "--generations"), type = "integer", default = 1, 
                help = "number of generations, for realistic local kinship", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_ind <- opt$n_ind
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
generations <- opt$generations

# name of directory containing data
name <- paste0(
    'sim',
    '-n', n_ind,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-g', generations
)

# move to where the data will be
# let's not overwrite things, under the assumption that the simulations are very expensive to redo
# if the output directory exists, assume all the files we want are there too.  Only a non-existent output directory will work
setwd( '../data/' )
if ( dir.exists( name ) )
    stop('Output exists, will not overwrite: ', name)
# else create directory and move into there
dir.create( name )
setwd( name )


############################
### POPULATION STRUCTURE ###
############################

# these are steps involving admixture and family structure abstractly, without genotypes at all (or other allele frequencies)

# define population structure
inbr_subpops <- 1 : k_subpops # FST values for K subpopulations
if (verbose)
    message('admix_prop_1d_linear')
obj <- admix_prop_1d_linear(
    n_ind = n_ind,
    k_subpops = k_subpops,
    bias_coeff = bias_coeff,
    coanc_subpops = inbr_subpops,
    fst = fst
)
# in this case return value is a named list with three items:
admix_proportions <- obj$admix_proportions # admixture proportions
inbr_subpops <- obj$coanc_subpops # rescaled Fst vector for intermediate subpops

# save bnpsd data to an RData file
save( admix_proportions, inbr_subpops, file = 'bnpsd.RData' )

if ( generations > 1 ) {
    # simulate realistic generations of families

    if (verbose)
        message('sim_children_generations_kinship')
    # defines parents semi-randomly based on avoidance of close relatives
    # but otherwise with strong assortative mating to preserve population structure
    data_G <- sim_children_generations_kinship(
        G = generations,
        n = n_ind,
        iterations = iterations,
        verbose = verbose
    )
    parents <- data_G$parents # list with a matrix per generation

    # save parents structure (for current code, it's most natural to use an RData file, rather than a FAM file, which we can write but can't parse back easily into the structure we want)
    save( parents, file = 'parents.RData' )
}
