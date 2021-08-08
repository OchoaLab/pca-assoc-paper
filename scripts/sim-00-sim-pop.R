# simulates population structure (without genotypes)
# based closely on ~/scripts/sim_reps_00_sim_pop.R on main scripts dir of this repo

# things that can't be changed:
# - Fst of intermediate subpopulations is a ramp proportional to 1:k
# - Admixture proportions are from 1D geography model

# creates:
# - bnpsd.RData

library(optparse) # for terminal options
library(bnpsd)    # simulate admixed population (structured population)
library(readr)    # to save FST value, for tables on paper

# hardcoded params
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

# save FST value of simulation
write_lines( fst, 'fst.txt' )

